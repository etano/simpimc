#ifndef BisectClass_H
#define BisectClass_H

#include "MoveClass.h"

class Bisect : public Move
{
private:
  string species;
  int iSpecies;
  int nImages;
  unsigned int nLevel, nBisectBeads;
  unsigned int bead0, bead1;
  double lambda;

  std::vector<Bead*> affBeads;
protected:

public:
  // Constructor
  Bisect(Path &tmpPath, RNG &tmpRNG, vector<Action*> &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual int Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif
