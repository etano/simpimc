#ifndef SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_

#include "observable_class.h"

class PairCorrelation : public Observable
{
private:
  std::string species_a, species_b;
  uint32_t species_a_i, species_b_i;
  Histogram gr;
protected:
public:
  PairCorrelation(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
    std::string data_type = "histogram";
    out.Write(prefix+"/data_type",data_type);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
