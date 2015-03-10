#ifndef VaryOptimizedNodalClass_H
#define VaryOptimizedNodalClass_H

#include "MoveClass.h"

class VaryOptimizedNodal : public Move
{
private:
  uint iSpecies, paramSet0, paramSet1;
protected:

public:
  // Constructor
  VaryOptimizedNodal(Path &tmpPath, RNG &tmpRNG, std::vector<std::shared_ptr<Action> > &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif
