#ifndef CoulombClass_H
#define CoulombClass_H

#include "ActionClass.h"

class Coulomb : public Action
{
private:
  int nImages;
  int nOrder;
  RealType Z1Z2;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int offsetA, offsetB;
  int maxLevel;
protected:

public:
  // Constructor
  Coulomb(Path &path, Input &in, IOClass &out)
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
