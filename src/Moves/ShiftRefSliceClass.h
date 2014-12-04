#ifndef ShiftRefSliceClass_H
#define ShiftRefSliceClass_H

#include "MoveClass.h"

class ShiftRefSlice : public Move
{
private:
  string species;
  int iSpecies;
  int refBead0, refBead1;
protected:

public:
  // Constructor
  ShiftRefSlice(Path &tmpPath, RNG &tmpRNG, vector<Action*> &actionList, Input &in, IOClass &out)
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
