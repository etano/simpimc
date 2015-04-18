#include "optimized_nodal_class.h"

// Initialize parameters
void OptimizedNodal::Init(Input &in)
{
  // Read in variational parameters
  model_i = in.GetAttribute<uint32_t>("model");
  std::vector<Input> param_set_inputs = in.GetChildList("ParamSet");
  for (auto& param_setInput : param_set_inputs) {
    std::vector<Input> param_inputs = param_setInput.GetChildList("Param");
    std::vector<double> params;
    for (auto& paramInput : param_inputs)
      params.push_back(paramInput.GetAttribute<double>("val"));
    param_sets.push_back(params);
  }

  // Read in common nodal things
  Nodal::Init(in);
}

// Create spline
void OptimizedNodal::SetupSpline()
{
  // Resize spline field
  uint32_t nSpline = path.n_bead/2 + (path.n_bead%2) + 1;

  // Create splines
  rho_free_r_splines.resize(param_sets.size());
  for (param_set_i=0; param_set_i<param_sets.size(); ++param_set_i) {
    for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
      double t_i_4_lambda_tau = Geti4LambdaTau(spline_i+1); // TODO: This is hard-coded for free-particle-like nodal structures.
      FreeSpline rho_free_r_spline(path.L, n_images, path.n_d, t_i_4_lambda_tau);
      rho_free_r_splines[param_set_i].push_back(rho_free_r_spline);
    }
    std::cout << "...param set " << param_set_i << " complete." << std::endl;
  }

  // Reset param_set_i
  SetParamSet(0);
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double OptimizedNodal::GetGij(const vec<double>& r, const uint32_t slice_diff)
{
  return rho_free_r_splines[param_set_i][slice_diff-1].GetGij(r);
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double OptimizedNodal::GetGijDGijDr(const vec<double>& r, const uint32_t slice_diff, vec<double>& dgij_dr)
{
  return rho_free_r_splines[param_set_i][slice_diff-1].GetGijDGijDr(r, dgij_dr);
}
