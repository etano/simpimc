#ifndef PermBisectClass_H
#define PermBisectClass_H

#include "MoveClass.h"

class PermBisect : public Move
{
private:
  unsigned int nLevel;
  unsigned int nBisectBeads;
  Tvector permTable;

  int DoPermBisect();
  RealType constructPermTable(const int bead0, const int bead1, const int nBisectBeads, const bool rollOver);
  int selectPerm(Ivector& permParts, RealType permTot);
  unsigned int permuteb(field<*Bead> b, int permType);

  std::vector<Bead*> affBeads;
protected:

public:
  PermBisect(Path& pathIn, RNG& rngIn, RealType perAcceptDesiredIn, int nEqSweepIn, int nEqStepIn, int moveSkipIn)
    : Move(pathIn, rngIn, perAcceptDesiredIn, nEqSweepIn, nEqStepIn, moveSkipIn )
  {
    moveLabel = "PermBisect";
    stepSize = floor(log(path.nBead/2.0)/log(2));

    // Initiate permutation table
    permTable.zeros(path.nPermType * path.nPart * (path.nPart-1) * (path.nPart-2));
  }

  virtual void MakeMove();
};

#endif
