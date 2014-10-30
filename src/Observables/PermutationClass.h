#ifndef PermutationClass_H
#define PermutationClass_H

#include "ObservableClass.h"

class Permutation : public Observable
{
private:
  Tvector cycles;
  vector<int> sectors;
  string species;
  int iSpecies, offset;
  bool firstSector;
protected:
public:
  Permutation(Path &tmpPath, Input &in, IOClass &out)
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
