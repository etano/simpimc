#ifndef SIMPIMC_ACTIONS_OPTIMIZED_FREE_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_OPTIMIZED_FREE_NODAL_CLASS_H_

#include "nodal_class.h"
#include "../../free_spline_class.h"

/// A nodal action class whose parameters may be varied
class OptimizedFreeNodal : public Nodal
{
private:
  uint32_t init_param_set; ///< Initial parameter set
  std::vector<std::vector<FreeSpline>> rho_free_splines; ///< Holds the splined action for every time slice and parameter set

  /// Returns the value of g_ij
  virtual double GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff)
  {
    return rho_free_splines[param_set_i][slice_diff-1].GetRhoFree(path.Dr(b_i,b_j));
  }

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double> &dgij_dr)
  {
    return rho_free_splines[param_set_i][slice_diff-1].GetGradRhoFree(path.Dr(b_i,b_j), dgij_dr);
  }

  // Returns 1/(4\lambda\tau)
  inline virtual double GetTau(const uint32_t slice_diff)
  {
    // Choose model
    switch(model_i) {
      case 0: // \alpha n \tau
        return path.tau*param_sets[param_set_i][0]*slice_diff;
        break;
      case 1: // \alpha[i] n \tau
        return path.tau*param_sets[param_set_i][slice_diff]*slice_diff;
      default: // n \tau
        return path.tau*slice_diff;;
        break;
    }
  };

  /// Creates splined action for all time slices and parameter sets
  virtual void SetupSpline()
  {
    // Create splines
    uint32_t n_spline = path.n_bead/2 + (path.n_bead%2) + 1;
    rho_free_splines.resize(param_sets.size());
    for (param_set_i=0; param_set_i<param_sets.size(); ++param_set_i) {
      for (uint32_t spline_i=0; spline_i<n_spline; ++spline_i)
        rho_free_splines[param_set_i].emplace_back(path.L, n_images, lambda, GetTau(spline_i+1), false);
      std::cout << "...param set " << param_set_i << " complete." << std::endl;
    }

    // Reset param_set_i
    SetParamSet(init_param_set);
  }

public:
  /// Constructor calls Init function
  OptimizedFreeNodal(Path &path, Input &in, IO &out)
    : Nodal(path,in,out)
  {
    // Read in variational parameters
    model_i = in.GetAttribute<uint32_t>("model");
    init_param_set = in.GetAttribute<uint32_t>("init_param_set");
    std::vector<Input> param_set_inputs = in.GetChildList("ParamSet");
    for (auto& param_setInput : param_set_inputs) {
      std::vector<Input> param_inputs = param_setInput.GetChildList("Param");
      std::vector<double> params;
      for (auto& paramInput : param_inputs)
        params.push_back(paramInput.GetAttribute<double>("val"));
      param_sets.push_back(params);
    }
    param_set_i = in.GetAttribute<uint32_t>("init_param_set",0);

    // Setup splines
    SetupSpline();

    // Test
    bool init_good = TestNodes();
    out.Write("Actions/"+name+"/init_good", init_good);
  }

  /// Write information about the action
  virtual void Write() {};

};

#endif // SIMPIMC_ACTIONS_OPTIMIZED_FREE_NODAL_CLASS_H_
