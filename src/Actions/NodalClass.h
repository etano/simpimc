#ifndef NodalClass_H
#define NodalClass_H

#include "ActionClass.h"

class Nodal : public Action
{
private:

protected:
  int nImages, maxLevel;
  string species;
  int iSpecies, offset;

public:
  // Constructor
  Nodal(Path &path, Input &in, IOClass &out)
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
