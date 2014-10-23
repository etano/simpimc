#ifndef PairCorrelationClass_H
#define PairCorrelationClass_H

#include "ObservableClass.h"

class PairCorrelation : public Observable
{
private:
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int offsetA, offsetB;
  Histogram gr;
protected:
public:
  PairCorrelation(Path &tmpPath, Input &in, IOClass &out)
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
