#ifndef SIMPIMC_ACTIONS_OPTIMIZED_SHO_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_OPTIMIZED_SHO_NODAL_CLASS_H_

#include "nodal_class.h"

/// A nodal action class whose parameters may be varied
class OptimizedSHONodal : public Nodal
{
private:
  /// Returns the value of g_ij
  virtual double GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff)
  {
    double omega = param_sets[param_set_i][0];
    double level_tau_omega = slice_diff*path.tau*omega;
    double cofactor = omega*i_4_lambda_tau/(M_PI*sinh(level_tau_omega));
    vec<double> r_i = path.GetR(b_i);
    vec<double> r_j = path.GetR(b_j);
    return pow(cofactor,path.n_d/2.)*exp(-cofactor*((dot(r_i,r_i) + dot(r_j,r_j))*cosh(level_tau_omega) - 2*dot(r_i,r_j)));
  }

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double> &dgij_dr)
  {
    double omega = param_sets[param_set_i][0];
    double level_tau_omega = slice_diff*path.tau*omega;
    vec<double> r_i = path.GetR(b_i);
    vec<double> r_j = path.GetR(b_j);
    double gij = GetGij(b_i, b_j, slice_diff);
    dgij_dr = -omega*i_4_lambda_tau/(M_PI*sinh(level_tau_omega))*(2.*r_i*cosh(level_tau_omega) - 2.*r_j)*gij;
    return gij;
  }

  /// Creates splined action for all time slices and parameter sets
  virtual void SetupSpline() {};
public:
  /// Constructor calls Init function
  OptimizedSHONodal(Path &path, Input &in, IO &out)
    : Nodal(path,in,out)
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
    param_set_i = in.GetAttribute<uint32_t>("init_param_set",0);

    // Test
    bool init_good = TestNodes();
    out.Write("Actions/"+name+"/init_good", init_good);
  }

  /// Write information about the action
  virtual void Write() {};

};

#endif // SIMPIMC_ACTIONS_OPTIMIZED_SHO_NODAL_CLASS_H_
