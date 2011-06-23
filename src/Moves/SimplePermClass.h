#ifndef SimplePermClass_H
#define SimplePermClass_H

#include "MoveClass.h"

class SimplePerm : public Move
{
private: 
  int DoSimplePerm( const int iBead );

  bool checkPermRadius( const int i , const int j , const int k );
  void setPermRadius( const int iBead );
  int permute( const int slice , const int i , const int j , const int k , int permType );
protected:

public:
  SimplePerm( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "Simple Perm";
    stepSize = 1.0;
  }

  virtual void MakeMove();
};

#endif
