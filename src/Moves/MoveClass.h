#ifndef MoveClass_H
#define MoveClass_H

#include "../StandardLibraries.h"       // Standard libraries
#include "../GlobalConstants.h"
#include "../PathClass.h"
#include "../BeadClass.h"
#include "../RNGClass.h"

class Move
{
private: 

protected:
  // Path
  Path& path;
  vec dr;

  // RNG
  RNG& rng;

  // Permutation Functions
  void assignParticleLabels();
  void setPerm( const int permType , int* perm , int* iPerm , const int i , const int j , const int k );

  // Equilibration
  double perAcceptDesired;
  unsigned int nEqSweep;
  unsigned int nEqStep;
  double stepSize;

  // Acceptance
  unsigned int nAttempt;
  unsigned int nAccept;
  double perAccept;
public:
  // Constructor
  Move( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : path(pathIn) , rng(rngIn) , perAcceptDesired(perAcceptDesiredIn) , nEqSweep(nEqSweepIn) , nEqStep(nEqStepIn) , moveSkip(moveSkipIn) 
  {
    dr.set_size(path.nD);
  }

  // Label
  const char* moveLabel;

  // Moves
  unsigned int moveSkip;
  virtual void MakeMove() {};

  // Equilibration
  void Equilibrate();

  // Acceptance
  void resetCounters();
  double getPerAccept();
};

#endif
