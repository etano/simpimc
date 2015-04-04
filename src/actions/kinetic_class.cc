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
  rho_free_r_splines.set_size(nSpline);

  // Create splines
  #pragma omp parallel for
  for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
    double t_i_4_lambda_tau = i_4_lambda_tau/(spline_i+1);

    // Make rho_free
    vec<double> rho_free_r, num_sum_r;
    rho_free_r.zeros(r_grid.num);
    if (spline_i == 0)
      num_sum_r.zeros(r_grid.num);
    for (uint32_t i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r2 = r*r;
      double r2_i_4_lambda_tau = r2*t_i_4_lambda_tau;
      for (uint32_t image=1; image<=n_images; image++) {
          double t_r_p(r+image*path.L);
          double exp_part_p = path.FastExp(r2_i_4_lambda_tau - t_r_p*t_r_p*t_i_4_lambda_tau);
          double t_r_m(r-image*path.L);
          double exp_part_m = path.FastExp(r2_i_4_lambda_tau - t_r_m*t_r_m*t_i_4_lambda_tau);
          rho_free_r(i) += exp_part_p + exp_part_m;
          if (spline_i == 0 && r2 != 0.)
            num_sum_r(i) += (t_r_p*t_r_p*exp_part_p + t_r_m*t_r_m*exp_part_m)/r2;
      }
      rho_free_r(i) = log1p(std::min(10.,rho_free_r(i)));
      if (spline_i == 0)
        num_sum_r(i) = log1p(std::min(10.,num_sum_r(i)));
    }
    BCtype_d xBC = {NATURAL, NATURAL};
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines(spline_i) = rho_free_r_spline;
    if (spline_i == 0)
      num_sum_r_spline = create_UBspline_1d_d(r_grid, xBC, num_sum_r.memptr());
  }

}

double Kinetic::GetGaussSum(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff)
{
  double gauss_sum;
  eval_UBspline_1d_d(rho_free_r_splines(slice_diff-1),r,&gauss_sum);
  return exp(-(r2_i_4_lambda_tau/slice_diff))*exp(0.9999*gauss_sum);
}

double Kinetic::GetGaussSumFast(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff)
{
  double gauss_sum;
  eval_UBspline_1d_d(rho_free_r_splines(slice_diff-1),r,&gauss_sum);
  return gauss_sum-(r2_i_4_lambda_tau/slice_diff);
};

double Kinetic::GetNumSum(const double r, const double r2_i_4_lambda_tau)
{
  double num_sum;
  eval_UBspline_1d_d(num_sum_r_spline,r,&num_sum);
  return -(r2_i_4_lambda_tau/path.tau)*exp(-r2_i_4_lambda_tau)*exp(0.9999*num_sum);
}

double Kinetic::DActionDBeta()
{
  double tot = n_part*path.n_bead*path.n_d/(2.*path.tau); // Constant term
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t p_i=0; p_i<n_part; p_i++) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      vec<double> num_sum(path.n_d), gauss_sum(path.n_d);
      vec<double> dr(path.Dr(path(species_i,p_i,b_i),path.GetNextBead(path(species_i,p_i,b_i),1)));
      double gauss_prod = 1.;
      for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
        double r2_i_4_lambda_tau = dr(d_i)*dr(d_i)*i_4_lambda_tau;
        num_sum(d_i) = GetNumSum(dr(d_i),r2_i_4_lambda_tau);
        gauss_sum(d_i) = GetGaussSum(dr(d_i),r2_i_4_lambda_tau,1);
        gauss_prod *= gauss_sum(d_i);
      }
      double scalar_num_sum = 0.;
      for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
        double num_prod = 1.;
        for (uint32_t d_j=0; d_j<path.n_d; d_j++) {
          if (d_i != d_j)
            num_prod *= gauss_sum(d_j);
          else
            num_prod *= num_sum(d_j);
        }
        scalar_num_sum += num_prod;
      }
      tot += scalar_num_sum/gauss_prod;
    }
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
        double gauss_prod_exp = 0;
        for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
          double r2_i_4_lambda_tau = dr(d_i)*dr(d_i)*i_4_lambda_tau;
          gauss_prod_exp += GetGaussSumFast(dr(d_i),r2_i_4_lambda_tau,skip);
        }
        tot -= gauss_prod_exp;
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
      double gauss_prod, gauss_sum, dist;
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
      double gauss_prod, gauss_sum, dist;
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
