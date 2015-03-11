#ifndef SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_

#include "observable_class.h"

class Permutation : public Observable
{
private:
  vec<double> cycles;
  std::vector<uint> sectors;
  std::string species;
  uint species_i;
  bool first_sector;
protected:
public:
  Permutation(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
