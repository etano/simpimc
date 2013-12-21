#ifndef EnergyClass_H
#define EnergyClass_H

#include "ObservableClass.h"

class Energy : public Observable
{
private:
  RealType E, KE, VE, NE;

  // Energy Observables
  RealType getKE();
  RealType getVE();
  RealType getVEint();
  RealType getVEext();
  RealType getNE();
protected:
public:
  Energy(Path &tmpPath, Input &in, IOClass &out)
    : Observable(tmpPath, in, out)
  {
    Init();
  }

  virtual void Init();
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
