#ifndef BisectClass_H
#define BisectClass_H

#include "MoveClass.h"

class Bisect : public Move
{
private:
  int DoBisect(const int iPart);
  std::vector<Bead*> affBeads;
protected:

public:
  // Constructor
  Bisect(Path &tmpPath, RNG &tmpRNG, vector<Action*> &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  int maxLevel;

  virtual void Init(Input &in);
  virtual void MakeMove();
};

#endif
