#include "free_spline_class.h"

// Create a spline for each possible slice_diff
void FreeSpline::Init(const double L, const uint32_t n_images)
{
  // Setup grid
  Ugrid r_grid;
  r_grid.start = -L/2.;
  r_grid.end = L/2.;
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Make rho_free_r_spline
  vec<double> rho_free_r;
  rho_free_r.zeros(r_grid.num);
  for (uint32_t i=0; i<r_grid.num; ++i) {
    double r = r_grid.start + i*dr;
    double r2_i_4_lambda_tau = r*r*i_4_lambda_tau;
    for (uint32_t image=1; image<=n_images; image++) {
      double t_r = r + image*L;
      rho_free_r(i) += exp(r2_i_4_lambda_tau - t_r*t_r*i_4_lambda_tau);
      t_r = r - image*L;
      rho_free_r(i) += exp(r2_i_4_lambda_tau - t_r*t_r*i_4_lambda_tau);
    }
    rho_free_r(i) = log1p(std::min(10.,rho_free_r(i)));
  }
  BCtype_d xBC = {NATURAL, NATURAL};
  rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeSpline::GetGij(const vec<double>& r)
{
  double gij = 1.;
  for (uint32_t d_i=0; d_i<n_d; d_i++) {
    double gij_image_action;
    eval_UBspline_1d_d(rho_free_r_spline,r(d_i),&gij_image_action);
    gij *= exp(gij_image_action - r(d_i)*r(d_i)*i_4_lambda_tau);
  }
  return gij;
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeSpline::GetGijDGijDr(const vec<double>& r, vec<double>& dgij_dr)
{
  double gij = 1.;
  for (uint32_t d_i=0; d_i<n_d; d_i++) {
    double gij_image_action, dgij_dr_image_action;
    eval_UBspline_1d_d_vg(rho_free_r_spline,r(d_i),&gij_image_action,&dgij_dr_image_action);
    double gij_d_i = exp(gij_image_action - r(d_i)*r(d_i)*i_4_lambda_tau);
    gij *= gij_d_i;
    dgij_dr(d_i) = dgij_dr_image_action - 2.*r(d_i)*i_4_lambda_tau;
  }
  dgij_dr *= gij;
  return gij;
}
