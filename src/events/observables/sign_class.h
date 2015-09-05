#ifndef SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
#define SIMPIMC_OBSERVABLES_SIGN_CLASS_H_

#include "observable_class.h"

/// Measure the sign weight
class Sign : public Observable
{
private:
  double sign; ///< Sign accumulator

  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    sign += path.CalcSign();
    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_measure = 0;
    sign = 0.;
  }
public:
  /// Constructor calls Init
  Sign(Path &path, Input &in, IO &out)
    : Observable(path, in, out, "scalar")
  {
    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      double norm = 1.*n_measure;

      // Write sign
      sign /= norm;
      if (first_time) {
        first_time = 0;
        out.CreateExtendableDataSet(prefix, "x", sign);
      } else {
        out.AppendDataSet(prefix, "x", sign);
      }

      Reset();
    }
  }

};

#endif // SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
