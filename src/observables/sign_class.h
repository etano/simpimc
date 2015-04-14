#ifndef SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
#define SIMPIMC_OBSERVABLES_SIGN_CLASS_H_

#include "observable_class.h"

/// Measure the sign weight
class Sign : public Observable
{
private:
  double sign; ///< Sign accumulator

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  Sign(Path &path, Input &in, IO &out)
    : Observable(path, in, out, "scalar")
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
