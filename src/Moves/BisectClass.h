#ifndef BisectClass_H
#define BisectClass_H

#include "MoveClass.h"

class Bisect : public Move
{
private: 
  int DoBisect( const int iPart );
  std::vector<Bead*> affBeads;
protected:

public:
  // Constructor
  Bisect( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Bisect";
    stepSize = 3.0;
  }

  virtual void MakeMove();
};

#endif
