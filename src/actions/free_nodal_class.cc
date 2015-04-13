#include "free_nodal_class.h"

// Create a spline for each possible slice_diff
void FreeNodal::SetupSpline()
{
  // Resize spline field
  uint32_t nSpline = path.n_bead/2 + (path.n_bead%2) + 1;

  // Create splines
  #pragma omp parallel for
  for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
    double t_i_4_lambda_tau = i_4_lambda_tau/(spline_i+1);
    FreeSpline rho_free_r_spline(path.L, n_images, path.n_d, t_i_4_lambda_tau);
    rho_free_r_splines.push_back(rho_free_r_spline);
  }
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGij(const vec<double>& r, const uint32_t slice_diff)
{
  return rho_free_r_splines[slice_diff-1].GetGij(r);
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGijDGijDr(const vec<double>& r, const uint32_t slice_diff, vec<double>& dgij_dr)
{
  return rho_free_r_splines[slice_diff-1].GetGijDGijDr(r, dgij_dr);
}
