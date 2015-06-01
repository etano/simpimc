#include "free_nodal_class.h"

// Create a spline for each possible slice_diff
void FreeNodal::SetupSpline()
{
  // Create splines
  uint32_t n_spline = path.n_bead/2 + (path.n_bead%2) + 1;
  #pragma omp parallel for
  for (uint32_t spline_i=0; spline_i<n_spline; ++spline_i) {
    rho_free_splines.emplace_back(path.L, n_images, lambda, path.tau*(spline_i+1), false);
  }
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff)
{
  return rho_free_splines[slice_diff-1].GetRhoFree(path.Dr(b_i,b_j));
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double FreeNodal::GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double>& dgij_dr)
{
  return rho_free_splines[slice_diff-1].GetGradRhoFree(path.Dr(b_i,b_j), dgij_dr);
}
