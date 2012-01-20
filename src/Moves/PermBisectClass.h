#ifndef PermBisectClass_H
#define PermBisectClass_H

#include "MoveClass.h"

class PermBisect : public Move
{
private: 
  unsigned int nLevel;
  unsigned int nBisectBeads;
  vec permTable;

  int DoPermBisect();
  double constructPermTable( const int bead0 , const int bead1 , const int nBisectBeads , const bool rollOver );
  int selectPerm( int* permParts , double permTot );
  unsigned int permuteb( Bead* b[3] , int permType );  
  
  std::vector<Bead*> affBeads;
protected:
   
public:
  PermBisect( Path& pathIn , RNG& rngIn , double perAcceptDesiredIn , int nEqSweepIn , int nEqStepIn , int moveSkipIn )
    : Move( pathIn , rngIn , perAcceptDesiredIn , nEqSweepIn , nEqStepIn , moveSkipIn )
  { 
    moveLabel = "PermBisect";
    stepSize = floor(log(path.nBead/2.0)/log(2));
    std::cerr << moveLabel << ": " << stepSize << endl;
    // Initiate permutation table
    permTable.zeros( path.nPermType * path.nPart * (path.nPart-1) * (path.nPart-2) );
  }

  virtual void MakeMove();
};

#endif
