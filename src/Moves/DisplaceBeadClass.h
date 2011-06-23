#ifndef DisplaceBeadClass_H
#define DisplaceBeadClass_H

#include "MoveClass.h"

class DisplaceBead : public Move
{
private: 
  int DoDisplaceBead( const int iPart , const int iBead );
  //int DoDisplace( const int iPart , const int iBead );
protected:
  
public:
  DisplaceBead( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Displace Bead";
    stepSize = 1.0;
  }

  virtual void MakeMove();
};

#endif
