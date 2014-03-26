#ifndef PairActionClass_H
#define PairActionClass_H

#include "ActionClass.h"

class PairAction : public Action
{
private:
  int nImages;
  int nOrder;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int offsetA, offsetB;
  int maxLevel;
protected:

public:
  // Constructor
  PairAction(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level);
  virtual void Write();
};

#endif
