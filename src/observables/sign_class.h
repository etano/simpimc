#ifndef SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
#define SIMPIMC_OBSERVABLES_SIGN_CLASS_H_

#include "observable_class.h"

class Sign : public Observable
{
private:
  double sign;
protected:
public:
  Sign(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_SIGN_CLASS_H_
