#ifndef SIMPIMC_ACTIONS_TRAP_CLASS_H_
#define SIMPIMC_ACTIONS_TRAP_CLASS_H_

#include "single_action_class.h"

/// Harmonic well action class
class Trap : public SingleAction
{
private:
  double omega; ///< Harmonic frequency of the trap
public:
  /// Constructor calls Init
  Trap(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    Init(in);
  }

  /// Initialize the action
  virtual void Init(Input &in);

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Write information about the action
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_TRAP_CLASS_H_
