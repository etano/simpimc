#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_

#include "observable_class.h"

/// Measures the contact density between two species of particles
///
/// Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
class ContactDensity : public Observable
{
private:
  double total; ///< Running total
  uint32_t z_a; ///< Charge of ion-like particle
  uint32_t species_a_i; ///< Ion species index
  uint32_t species_b_i; ///< Other species index
  std::string species_a; ///< Name of ion species
  std::string species_b; ///< Name of other species
  std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that involves species_a or species_b
  std::vector<std::shared_ptr<Action>> &full_action_list; ///< Vector of pointers to all actions

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  ContactDensity(Path &path, std::vector<std::shared_ptr<Action>>& t_action_list, Input &in, IO &out)
    : full_action_list(t_action_list), Observable(path, in, out, "scalar")
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
