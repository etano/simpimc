#include "free_nodal_class.h"

// Create a spline for each possible slice_diff
void FreeNodal::SetupSpline()
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
    vec<double> rho_free_r;
    rho_free_r.zeros(r_grid.num);
    for (uint32_t i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r2_i_4_lambda_tau = r*r*t_i_4_lambda_tau;
      for (uint32_t image=1; image<=n_images; image++) {
        double t_r = r + image*path.L;
        rho_free_r(i) += path.FastExp(r2_i_4_lambda_tau - t_r*t_r*t_i_4_lambda_tau);
        t_r = r - image*path.L;
        rho_free_r(i) += path.FastExp(r2_i_4_lambda_tau - t_r*t_r*t_i_4_lambda_tau);
      }
      rho_free_r(i) = log1p(std::min(10.,rho_free_r(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines(spline_i) = rho_free_r_spline;
  }
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGij(const vec<double>& r, const uint32_t slice_diff)
{
  double gij = 1.;
  double i_4_lambda_level_tau = i_4_lambda_tau/slice_diff;
  for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
    double gij_image_action;
    eval_UBspline_1d_d(rho_free_r_splines(slice_diff-1),r(d_i),&gij_image_action);
    gij *= exp(0.9999*gij_image_action)*exp(-(r(d_i)*r(d_i)*i_4_lambda_level_tau));
  }
  return gij;
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGijDGijDr(const vec<double>& r, const uint32_t slice_diff, vec<double>& dgij_dr)
{
  double gij = 1.;
  double i_4_lambda_level_tau = i_4_lambda_tau/slice_diff;
  for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
    double gij_image_action, dgij_dr_image_action;
    eval_UBspline_1d_d_vg(rho_free_r_splines(slice_diff-1),r(d_i),&gij_image_action,&dgij_dr_image_action);
    double gij_d_i = exp(0.9999*gij_image_action)*exp(-(r(d_i)*r(d_i)*i_4_lambda_level_tau));
    gij *= gij_d_i;
    dgij_dr(d_i) = (dgij_dr_image_action - r(d_i)*i_4_lambda_level_tau)*gij_d_i;
  }
  return gij;
}
