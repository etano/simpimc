#include "free_nodal_class.h"

// Initialize parameters
void FreeNodal::Init(Input &in)
{
  // Read in things
  n_images = in.GetAttribute<int>("n_images");
  species = in.GetAttribute<std::string>("species");
  species_list.push_back(species);
  std::cout << "Setting up nodal action for " << species << "..." << std::endl;
  max_level = in.GetAttribute<uint32_t>("max_level",0);
  path.GetSpeciesInfo(species,species_i);
  n_part = path.species_list[species_i]->n_part;
  i_4_lambda_tau = 1./(4.*path.species_list[species_i]->lambda*path.tau);

  // Write things to file
  out.Write("Actions/"+name+"/n_images", n_images);
  out.Write("Actions/"+name+"/species", species);
  out.Write("Actions/"+name+"/max_level", max_level);

  // Setup splines
  SetupSpline();

  // Set up determinants
  rho_f.set_size(path.n_bead);
  rho_f_c.set_size(path.n_bead);

  // Test initial configuration
  std::vector< std::pair<uint32_t,uint32_t>> particles;
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    particles.push_back(std::make_pair(species_i,p_i));
  bool init_good = 1;
  path.SetMode(1);
  if (GetAction(0, path.n_bead, particles, 0) == 1.e100) {
    std::cout << "Warning: initializing with broken nodes!" << std::endl;
    init_good = 0;
  }
  out.Write("Actions/"+name+"/init_good", init_good);
}

// Create a spline for each possible slice_diff
void FreeNodal::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  if (path.pbc) {
    r_grid.start = -path.L/2.;
    r_grid.end = path.L/2.;
  } else {
    r_grid.start = 100.;
    r_grid.end = 100.;
    n_images = 0;
  }
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  uint32_t nSpline = path.n_bead/2 + (path.n_bead%2) + 1;
  rho_free_r_splines.set_size(nSpline);

  // Create splines
  for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
    vec<double> rho_free_r(r_grid.num);
    double t_i_4_lambda_tau = i_4_lambda_tau/(spline_i+1);

    // Make rho_free
    for (uint32_t i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r2_i_4_lambda_tau = r*r*t_i_4_lambda_tau;
      rho_free_r(i) = 0.;
      for (int image=-n_images; image<=n_images; image++) {
        if (image != 0) {
          double t_r = r + image*path.L;
          rho_free_r(i) += path.FastExp(r2_i_4_lambda_tau - t_r*t_r*t_i_4_lambda_tau);
        }
      }
      rho_free_r(i) = log1p(std::min(10.,rho_free_r(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines(spline_i) = rho_free_r_spline;
  }
}

// Evaluate \rho_{ij}
double FreeNodal::GetGij(const vec<double>& r, const uint32_t slice_diff)
{
  double gauss_prod = 1.;
  for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
    double gauss_sum;
    eval_UBspline_1d_d(rho_free_r_splines(slice_diff-1),r(d_i),&gauss_sum);
    gauss_sum = exp(0.9999*gauss_sum);
    gauss_sum *= exp(-(r(d_i)*r(d_i)*i_4_lambda_tau/slice_diff));
    gauss_prod *= gauss_sum;
  }
  return gauss_prod;
}
