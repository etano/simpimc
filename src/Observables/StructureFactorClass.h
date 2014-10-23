#ifndef StructureFactorClass_H
#define StructureFactorClass_H

#include "ObservableClass.h"

class StructureFactor : public Observable
{
private:
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int offsetA, offsetB;
  Tvector sk;
  RealType kCut;
protected:
public:
  StructureFactor(Path &tmpPath, Input &in, IOClass &out)
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
