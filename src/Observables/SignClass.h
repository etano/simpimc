#ifndef SignClass_H
#define SignClass_H

#include "ObservableClass.h"

class Sign : public Observable
{
private:
  double sign;
protected:
public:
  Sign(Path &tmpPath, Input &in, IOClass &out)
    : Observable(tmpPath, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
