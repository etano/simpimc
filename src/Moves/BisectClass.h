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
    stepSize = floor(log(path.nBead/2.0)/log(2));
  }

  // Make move
  virtual void MakeMove();
  virtual void Write();
};

#endif
