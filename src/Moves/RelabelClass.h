#ifndef RelabelClass_H
#define RelabelClass_H

#include "MoveClass.h"

class Relabel : public Move
{
private: 
  int DoRelabel();
protected:

public:
  Relabel( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Relabel";
    stepSize = 1.0;
  }

  virtual void MakeMove();
};

#endif
