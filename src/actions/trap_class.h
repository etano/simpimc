#ifndef SIMPIMC_ACTIONS_TRAP_CLASS_H_
#define SIMPIMC_ACTIONS_TRAP_CLASS_H_

#include "single_action_class.h"

/// Harmonic well action class
class Trap : public SingleAction
{
private:
  double omega; ///< Harmonic frequency of the trap
  double cofactor_a; ///< Constant from 4th order expansion of the action
  double cofactor_b; ///< Constant from 4th order expansion of the beta derivative of the action
public:
  /// Constructor calls Init
  Trap(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    omega = in.GetAttribute<double>("omega");
    out.Write("/Actions/"+name+"/omega", omega);

    cofactor_a = 0.5*path.tau*omega*omega*(1. + path.tau*path.tau*omega*omega/12.);
    cofactor_b = 0.5*omega*omega*(1. + 3.*path.tau*path.tau*omega*omega/12.);
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Write information about the action
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_TRAP_CLASS_H_
