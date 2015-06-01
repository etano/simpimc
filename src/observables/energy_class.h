#ifndef SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
#define SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_

#include "observable_class.h"

/// Measure the energy of each action
class Energy : public Observable
{
private:
  bool first_sector; ///< Whether or not this is the first sector being written to file
  bool measure_per_sector; ///< Whether or not to measure per permutation section
  bool measure_potential; ///< Whether or not to measure the potential
  bool use_thermal_estimator; ///< Whether or not to use the thermal estimator of the energy
  bool use_virial_estimator; ///< Whether or not to use the virial estimator of the energy
  uint32_t b_i_last; ///< Which bead to average over for virial estimator
  uint32_t virial_window_size; ///< Size of averaging window for virial estimator
  std::vector<std::vector<std::pair<uint32_t,double>>> sector_energies; ///< Vector of vectors of (sector,energy_per_sector) pairs for each species
  std::vector<std::shared_ptr<Action>> &action_list; ///< Reference to vector of all actions
  vec<double> energies; ///< Vector of energies for each action
  vec<double> potentials; ///< Vector of potentials for each action

  /// Thermal estimator of the energy
  double ThermalEstimator();

  /// Virial estimator of the energy
  double VirialEstimator();

  /// Potential estimator
  double PotentialEstimator();

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
