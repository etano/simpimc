#ifndef SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_
#define SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_

#include "observable_class.h"

/// Record the structure factor between two species of particles
class StructureFactor : public Observable
{
private:
  double k_cut; ///< Cutoff in k space
  uint32_t species_a_i; ///< Index of first species
  uint32_t species_b_i; ///< Index of second species
  std::string species_a; ///< Name of first species
  std::string species_b; ///< Name of second species
  vec<double> sk; ///< Vector representing s(k)

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init and sets output data_type
  StructureFactor(Path &path, Input &in, IO &out)
    : Observable(path, in, out, "histogram")
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_
