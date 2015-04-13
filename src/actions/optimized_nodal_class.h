#ifndef SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_

#include "nodal_class.h"

/// A nodal action class whose parameters may be varied
class OptimizedNodal : public Nodal
{
private:

protected:
  field<UBspline_1d_d*> rho_node_r_splines; ///< Holds the splined action for every time slice and parameter set

  /// Returns the value of g_ij
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff);

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const vec<double> &r, const uint32_t slice_diff, vec<double> &dgij_dr);

  // Returns 1/(4\lambda\tau)
  inline virtual double Geti4LambdaTau(const uint32_t slice_diff)
  {
    // Choose model
    switch(model_i) {
      case 0: // \alpha/(4\lambda n\tau)
        return i_4_lambda_tau*param_sets[param_set_i][0]/slice_diff;
        break;
      case 1: // \alpha[i]/(4\lambda n\tau)
        return i_4_lambda_tau*param_sets[param_set_i][slice_diff]/slice_diff;
      default: // 1/(4\lambda n\tau)
        return i_4_lambda_tau/slice_diff;
        break;
    }
  };

  /// Creates splined action for all time slices and parameter sets
  virtual void SetupSpline();

public:
  /// Constructor calls Init function
  OptimizedNodal(Path &path, Input &in, IO &out)
    : Nodal(path,in,out)
  {
    Init(in);
  }

  /// Initialize the action
  virtual void Init(Input &in);

  /// Write information about the action
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_
