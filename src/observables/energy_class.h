#ifndef SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
#define SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_

#include "observable_class.h"

/// Measure the energy of each action
class Energy : public Observable
{
private:
  bool measure_potential; ///< Whether or not to measure the potential
  bool measure_per_sector; ///< Whether or not to measure per permutation section
  bool first_sector; ///< Whether or not this is the first sector being written to file
  uint32_t species_i; ///< Species index for per permutation sector measuring
  std::vector<std::pair<uint32_t,double>> sector_energies; ///< Vector of (sector,energy_per_sector) pairs
  std::vector<std::shared_ptr<Action>> &action_list; ///< Reference to vector of all actions
  vec<double> energies; ///< Vector of energies for each action
  vec<double> potentials; ///< Vector of potentials for each action

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  Energy(Path &path, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &out)
    : action_list(tmp_action_list), Observable(path, in, out)
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
