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
  rho_node_r_splines.set_size(param_sets.size(), nSpline);

  // Create splines
  for (uint32_t param_set_i=0; param_set_i<param_sets.size(); ++param_set_i) {
    #pragma omp parallel for
    for (uint32_t spline_i=0; spline_i<nSpline; ++spline_i) {
      vec<double> rho_node_r(r_grid.num);
      double t_i_4_lambda_tau = Geti4LambdaTau(spline_i+1); // TODO: This is hard-coded for free-particle-like nodal structures.

      // Make rho_free
      for (uint32_t i=0; i<r_grid.num; ++i) {
        double r = r_grid.start + i*dr;
        double r2i_4_lambda_tau = r*r*t_i_4_lambda_tau;
        rho_node_r(i) = 0.;
        for (int image=-n_images; image<=n_images; image++) {
          if (image != 0) {
            double t_r = r + image*path.L;
            rho_node_r(i) += path.FastExp(r2i_4_lambda_tau - t_r*t_r*t_i_4_lambda_tau);
          }
        }
        rho_node_r(i) = log1p(std::min(10.,rho_node_r(i)));
      }
      BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
      UBspline_1d_d* rho_node_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_node_r.memptr());
      rho_node_r_splines(param_set_i,spline_i) = rho_node_r_spline;
    }
    std::cout << "...param set " << param_set_i << " complete." << std::endl;
  }

  // Reset param_set_i
  SetParamSet(0);
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double OptimizedNodal::GetGij(const vec<double>& r, const uint32_t slice_diff)
{
  double gij = 1.;
  double i_4_lambda_level_tau = Geti4LambdaTau(slice_diff);
  for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
    double gij_image_action;
    eval_UBspline_1d_d(rho_node_r_splines(param_set_i,slice_diff-1),r(d_i),&gij_image_action);
    gij *= exp(0.9999*gij_image_action)*exp(-(r(d_i)*r(d_i)*i_4_lambda_level_tau));
  }
  return gij;
}

// Evaluate \rho_{ij} and d\rho_{ij}/dr_{ij}
double OptimizedNodal::GetGijDGijDr(const vec<double>& r, const uint32_t slice_diff, vec<double>& dgij_dr)
{
  double gij = 1.;
  double i_4_lambda_level_tau = Geti4LambdaTau(slice_diff);
  for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
    double gij_image_action, dgij_dr_image_action;
    eval_UBspline_1d_d_vg(rho_node_r_splines(param_set_i,slice_diff-1),r(d_i),&gij_image_action,&dgij_dr_image_action);
    double gij_d_i = exp(0.9999*gij_image_action)*exp(-(r(d_i)*r(d_i)*i_4_lambda_level_tau));
    gij *= gij_d_i;
    dgij_dr(d_i) = (dgij_dr_image_action - r(d_i)*i_4_lambda_level_tau)*gij_d_i;
  }
  return gij;
}

double OptimizedNodal::Geti4LambdaTau(const uint32_t slice_diff)
{
  double t_i_4_lambda_tau(i_4_lambda_tau);

  // Choose model
  switch(model_i) {
    case 0:
      t_i_4_lambda_tau *= param_sets[param_set_i][0]/slice_diff;
      break;
    default:
      t_i_4_lambda_tau /= slice_diff;
      break;
  }

  return t_i_4_lambda_tau;
}

void OptimizedNodal::Write() {}
