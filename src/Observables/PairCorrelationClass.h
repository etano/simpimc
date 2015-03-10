#ifndef PairCorrelationClass_H
#define PairCorrelationClass_H

#include "ObservableClass.h"

class PairCorrelation : public Observable
{
private:
  string speciesA, speciesB;
  uint iSpeciesA, iSpeciesB;
  Histogram gr;
protected:
public:
  PairCorrelation(Path &tmpPath, Input &in, IOClass &out)
    : Observable(tmpPath, in, out)
  {
    Init(in);
    string data_type = "histogram";
    out.Write(prefix+"/data_type",data_type);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
