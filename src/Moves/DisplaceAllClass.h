#ifndef DisplaceAllClass_H
#define DisplaceAllClass_H

#include "MoveClass.h"

class DisplaceAll : public Move
{
private: 
  int DoDisplaceAll();
protected:

public:
  DisplaceAll( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Displace All";
    stepSize = 1.0;
  }

  virtual void MakeMove();
};

#endif
