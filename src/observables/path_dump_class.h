#ifndef SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
#define SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_

#include "observable_class.h"

/// Dump out all path information
class PathDump : public Observable
{
private:
  uint32_t n_dump; ///< Number of times dumped
  uint32_t n_write_calls; ///< Number of times write has been called

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  PathDump(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
