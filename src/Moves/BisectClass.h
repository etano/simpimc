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
    stepSize = floor(log(path.nBead/2.0)/log(2));
  }

  virtual void MakeMove();
};

#endif
