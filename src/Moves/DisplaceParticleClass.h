#ifndef DisplaceParticleClass_H
#define DisplaceParticleClass_H

#include "MoveClass.h"

class DisplaceParticle : public Move
{
private:
  string species;
  int iSpecies;
  vector< std::shared_ptr<Bead> > affBeads;
  double stepSize;
protected:

public:
  // Constructor
  DisplaceParticle(Path &tmpPath, RNG &tmpRNG, std::vector< std::shared_ptr<Action> > &actionList, Input &in, IOClass &out)
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
