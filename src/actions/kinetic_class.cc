#include "kinetic_class.h"

void Kinetic::Init(Input &in)
{
  // Read in things
  n_images = in.GetAttribute<int>("n_images");
  species = in.GetAttribute<std::string>("species");
  species_list.push_back(species);
  std::cout << "Setting up kinetic action for " << species << "..." << std::endl;
  path.GetSpeciesInfo(species,species_i);
  n_part = path.species_list[species_i]->n_part;
  i_4_lambda_tau = 1./(4.*path.species_list[species_i]->lambda*path.tau);

  // Write things to file
  out.Write("Actions/"+name+"/n_images", n_images);
  out.Write("Actions/"+name+"/species", species);

  // Setup spline
  SetupSpline();
}

// Create a spline for each possible slice_diff
// TODO: Combine with free nodal action splines
void Kinetic::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  if (path.pbc) {
    r_grid.start = -path.L/2.;
    r_grid.end = path.L/2.;
  } else {
    r_grid.start = -100.;
    r_grid.end = 100.;
    n_images = 0;
  }
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  uint32_t nSpline = path.n_bead/2 + (path.n_bead%2) + 1;
  image_action_splines.set_size(nSpline);

  // Create splines
  #pragma omp parallel for
  for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
    double t_i_4_lambda_tau = i_4_lambda_tau/(spline_i+1);

    // Make image_action, d_image_action_d_tau, d_image_action_d_r
    vec<double> image_action, d_image_action_d_tau, d_image_action_d_r;
    image_action.zeros(r_grid.num);
    if (spline_i == 0) {
      d_image_action_d_tau.zeros(r_grid.num);
      d_image_action_d_r.zeros(r_grid.num);
    }
    for (uint32_t i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r2 = r*r;
      double r2_i_4_lambda_tau = r2*t_i_4_lambda_tau;
      for (uint32_t image=1; image<=n_images; image++) {
        double r_p(r+image*path.L);
        double r_p_2_i_4_lambda_tau = r_p*r_p*t_i_4_lambda_tau;
        double d_r_2_r_p_2_i_4_lambda_tau = r2_i_4_lambda_tau - r_p_2_i_4_lambda_tau;
        double exp_part_p = path.FastExp(d_r_2_r_p_2_i_4_lambda_tau);
        double r_m(r-image*path.L);
        double r_m_2_i_4_lambda_tau = r_m*r_m*t_i_4_lambda_tau;
        double d_r_2_r_m_2_i_4_lambda_tau = r2_i_4_lambda_tau - r_m_2_i_4_lambda_tau;
        double exp_part_m = path.FastExp(d_r_2_r_m_2_i_4_lambda_tau);
        image_action(i) += exp_part_p + exp_part_m;
        if (spline_i == 0) {
          d_image_action_d_r(i) += ((r-r_p)*exp_part_p + (r-r_m)*exp_part_m)*2.*t_i_4_lambda_tau;
          d_image_action_d_tau(i) += (d_r_2_r_p_2_i_4_lambda_tau*exp_part_p + d_r_2_r_m_2_i_4_lambda_tau*exp_part_m)/path.tau;
        }
      }
      if (spline_i == 0) {
        d_image_action_d_r(i) = -d_image_action_d_r(i)/(1.+image_action(i));
        d_image_action_d_tau(i) = d_image_action_d_tau(i)/(1.+image_action(i));
        //std::cout << r << " " << -log1p(std::min(10.,image_action(i))) << " " << image_action(i) << " " << d_image_action_d_r(i) << " " << d_image_action_d_tau(i) << std::endl;
      }
      image_action(i) = -log1p(image_action(i));
    }
    BCtype_d xBC = {NATURAL, NATURAL};
    UBspline_1d_d* image_action_spline = create_UBspline_1d_d(r_grid, xBC, image_action.memptr());
    image_action_splines(spline_i) = image_action_spline;
    if (spline_i == 0) {
      d_image_action_d_tau_spline = create_UBspline_1d_d(r_grid, xBC, d_image_action_d_tau.memptr());
      d_image_action_d_r_spline = create_UBspline_1d_d(r_grid, xBC, d_image_action_d_r.memptr());
    }
  }

}

double Kinetic::GetRhoFree(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff)
{
  return exp(GetLogRhoFree(r,r2_i_4_lambda_tau,slice_diff));
}

double Kinetic::GetLogRhoFree(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff)
{
  double image_action;
  eval_UBspline_1d_d(image_action_splines(slice_diff-1),r,&image_action);
  return -(image_action + r2_i_4_lambda_tau/slice_diff);
};

double Kinetic::GetDRhoFreeDTau(const double r, const double r2_i_4_lambda_tau)
{
  return GetDLogRhoFreeDTau(r,r2_i_4_lambda_tau)*GetRhoFree(r,r2_i_4_lambda_tau,1);
}

double Kinetic::GetDLogRhoFreeDTau(const double r, const double r2_i_4_lambda_tau)
{
  double d_image_action_d_tau;
  eval_UBspline_1d_d(d_image_action_d_tau_spline,r,&d_image_action_d_tau);
  return -((r2_i_4_lambda_tau/path.tau) + d_image_action_d_tau);
}

double Kinetic::GetDRhoFreeDR(const double r, const double r2_i_4_lambda_tau)
{
  return GetDLogRhoFreeDR(r,r2_i_4_lambda_tau)*GetRhoFree(r,r2_i_4_lambda_tau,1);
}

double Kinetic::GetDLogRhoFreeDR(const double r, const double r2_i_4_lambda_tau)
{
  double d_image_action_d_r;
  eval_UBspline_1d_d(d_image_action_d_r_spline,r,&d_image_action_d_r);
  return -(r*.2*i_4_lambda_tau + d_image_action_d_r);
}

double Kinetic::DActionDBeta()
{
  double tot = n_part*path.n_bead*path.n_d/(2.*path.tau); // Constant term
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
    for (uint32_t p_i=0; p_i<n_part; p_i++) {
      vec<double> dr(path.Dr(path(species_i,p_i,b_i),path.GetNextBead(path(species_i,p_i,b_i),1)));
      for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
        double r2_i_4_lambda_tau = dr(d_i)*dr(d_i)*i_4_lambda_tau;
        tot += GetDLogRhoFreeDTau(dr(d_i),r2_i_4_lambda_tau);
      }
    }
  }

  return tot;
}

double Kinetic::VirialEnergy(const double virial_window_size)
{
  // Constant term
  double tot = n_part*path.n_bead*path.n_d/(2.*virial_window_size*path.tau);

  // Permutation/winding term
  for (uint32_t p_i=0; p_i<n_part; p_i++) {

    // First bead
    std::shared_ptr<Bead> b_0(path(species_i,p_i,0));
    std::shared_ptr<Bead> b_1(path.GetNextBead(b_0,1));
    std::shared_ptr<Bead> b_2(path.GetNextBead(b_1,1));

    vec<double> dr_i_1(path.Dr(b_1,b_0));
    vec<double> dr_i_L(dr_i_1);
    for (uint32_t skip=1; skip<virial_window_size; skip++) {
      dr_i_L += path.Dr(b_2, b_1);
      b_1 = b_2;
      b_2 = path.GetNextBead(b_2,1);
    }
    tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.tau);

    // Other beads
    for (uint32_t b_i=1; b_i<path.n_bead; b_i++) {
      // Subtract first link and add new last link
      dr_i_L -= dr_i_1;
      dr_i_L += path.Dr(b_2,b_1);

      // Move beads one forward
      b_0 = path.GetNextBead(b_0,1);
      dr_i_1 = path.Dr(path.GetNextBead(b_0,1),b_0);
      b_1 = b_2;
      b_2 = path.GetNextBead(b_2,1);

      // Add to total
      tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.tau);
    }

      //// Image tau derivative
      //vec<double> dr_i(path.Dr(b_0, path.GetPrevBead(b_0,1)));
      //for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
      //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
      //  tot -= GetDLogRhoFreeDTau(dr_i(d_i),r2_i_4_lambda_tau) + r2_i_4_lambda_tau/path.tau;
      //}

      //// Image r derivative
      //vec<double> image_action_gradient(path.n_d);
      //for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
      //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
      //  image_action_gradient(d_i) = GetDLogRhoFreeDR(dr_i(d_i),r2_i_4_lambda_tau) + dr_i(d_i)*.2*i_4_lambda_tau;
      //  r2_i_4_lambda_tau = dr_i_1(d_i)*dr_i_1(d_i)*i_4_lambda_tau;
      //  image_action_gradient(d_i) += GetDLogRhoFreeDR(dr_i_1(d_i),r2_i_4_lambda_tau) + dr_i_1(d_i)*.2*i_4_lambda_tau;
      //}
      //tot -= dot(image_action_gradient,path.Dr(b_0,path(species_i,p_i,0)))/(2.*path.tau);

  }
  return tot;
}

double Kinetic::GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  uint32_t skip = 1<<level;
  double i_4_lambda_level_tau = i_4_lambda_tau/skip;
  double tot = 0.;
  for (auto& p: particles) {
    uint32_t s_i = p.first;
    uint32_t p_i = p.second;
    if (s_i == species_i) {
      std::shared_ptr<Bead> beadA(path(s_i,p_i,b0));
      std::shared_ptr<Bead> beadF(path.GetNextBead(beadA,b1-b0));
      while(beadA != beadF) {
        std::shared_ptr<Bead> beadB(path.GetNextBead(beadA,skip));
        vec<double> dr(path.Dr(beadA,beadB));
        for (uint32_t d_i=0; d_i<path.n_d; d_i++)
          tot -= GetLogRhoFree(dr(d_i),dr(d_i)*dr(d_i)*i_4_lambda_tau,skip);
        beadA = beadB;
      }
    }
  }

  return tot;
}

vec<double> Kinetic::GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  uint32_t skip = 1<<level;
  double i_4_lambda_level_tau = i_4_lambda_tau/skip;
  vec<double> tot;
  tot.zeros(path.n_d);
  std::shared_ptr<Bead> beadA, beadB, beadC, beadF;
  for (auto& p: particles) {
    uint32_t s_i = p.first;
    uint32_t p_i = p.second;
    if (s_i == species_i) {
      double gauss_prod, rho_free, dist;
      beadA = path(species_i,p_i,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        beadB = path.GetPrevBead(beadA,skip);
        vec<double> dr(path.Dr(beadB,beadA));
        tot -= dr;
        beadC = path.GetNextBead(beadA,skip);
        dr = path.Dr(beadA,beadC);
        tot += dr;
        beadA = beadC;
      }
    }
  }

  return 2.*i_4_lambda_level_tau*tot;
}

double Kinetic::GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  uint32_t skip = 1<<level;
  double i_4_lambda_level_tau = i_4_lambda_tau/skip;
  double tot = 0.;
  vec<double> dr(path.n_d);
  std::shared_ptr<Bead> beadA, beadF;
  for (auto& p: particles) {
    uint32_t s_i = p.first;
    uint32_t p_i = p.second;
    if (s_i == species_i) {
      double gauss_prod, rho_free, dist;
      beadA = path(species_i,p_i,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        tot += path.n_d*4.*i_4_lambda_level_tau;
        beadA = path.GetNextBead(beadA,skip);
      }
    }
  }

  return tot;
}

void Kinetic::Write()
{

}
