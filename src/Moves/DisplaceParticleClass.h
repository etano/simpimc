#ifndef DisplaceParticleClass_H
#define DisplaceParticleClass_H

#include "MoveClass.h"

class DisplaceParticle : public Move
{
private: 
  int DoDisplaceParticle( const int iPart );
protected:

public:
  DisplaceParticle( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Displace Particle";
    stepSize = 1.0;
  }

  virtual void MakeMove();
};

#endif
