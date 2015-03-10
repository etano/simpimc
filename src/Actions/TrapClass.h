#ifndef TrapClass_H
#define TrapClass_H

#include "ActionClass.h"

class Trap : public Action
{
private:
  int nImages;
  uint iSpecies, maxLevel;
  double omega;
  string species;
protected:

public:
  // Constructor
  Trap(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const vector< pair<uint,uint> > &particles, const uint level);
  virtual void Write();
};

#endif
