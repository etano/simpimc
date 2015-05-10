#ifndef SIMPIMC_ACTIONS_KINETIC_CLASS_H_
#define SIMPIMC_ACTIONS_KINETIC_CLASS_H_

#include "single_action_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

/// A kinetic action class
class Kinetic : public SingleAction
{
private:
  field<UBspline_1d_d*> image_action_splines; ///< Holds the splined action for every time slice
  UBspline_1d_d* d_image_action_d_tau_spline; ///< Holds the num_sum used in the calculation of the kinetic energy
  UBspline_1d_d* d_image_action_d_r_spline; ///< Holds the num_sum used in the calculation of the kinetic energy

  /// Creates splined action for all time slices
  void SetupSpline();

  /// Returns the free particle propagator
  double GetRhoFree(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff);

  /// Returns log of the free particle propagator
  double GetLogRhoFree(const double r, const double r2_i_4_lambda_tau, const uint32_t slice_diff);

  /// Returns the tau derivative of the free particle propagator
  double GetDRhoFreeDTau(const double r, const double r2_i_4_lambda_tau);

  /// Returns the tau derivative of the log of the free particle propagator
  double GetDLogRhoFreeDTau(const double r, const double r2_i_4_lambda_tau);

  /// Returns the spatial derivative of the free particle propagator
  double GetDRhoFreeDR(const double r, const double r2_i_4_lambda_tau);

  /// Returns the spatial derivative of the log of the free particle propagator
  double GetDLogRhoFreeDR(const double r, const double r2_i_4_lambda_tau);
public:
  /// Constructor calls Init
  Kinetic(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    Init(in);
  }

  /// Initialize the action
  virtual void Init(Input &in);

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the virial contribution of the action
  virtual double VirialEnergy(const double virial_window_size);

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Write information about the action
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_KINETIC_CLASS_H_
