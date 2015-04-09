#ifndef SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_

#include "observable_class.h"

/// Tracks the permutation sectors and cycles
class Permutation : public Observable
{
private:
  bool first_sector; ///< Whether or not the first permutation sector has been written
  uint32_t species_i; ///< Index of relevant species
  std::string species; ///< Name of relevant species
  std::vector<uint32_t> sectors; ///< Vector of sector indices
  vec<double> cycles; ///< Vector of cycle counts

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  Permutation(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
