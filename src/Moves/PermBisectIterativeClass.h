#ifndef PermBisectIterativeClass_H
#define PermBisectIterativeClass_H

#include "MoveClass.h"

class PermBisectIterative : public Move
{
private:
  string species;
  int iSpecies;
  int offset;
  int nImages;
  unsigned int nLevel, nBisectBeads, nPart, nPermPart, nPermType;
  RealType lambda, i4LambdaTauNBisectBeads, epsilon, logEpsilon;
  bool fixedNode;

  struct Cycle
  {
    RealType weight, contribution;
    int index, type;
    Ivector perm, iPerm, part;
  };
  vector<Cycle*> cycles;
  field<Cycle> all_cycles;
  Tmatrix t;

  Ivector permAttempt, permAccept;

  int DoPermBisect();
  void updatePermTable(const int bead0, const int nBisectBeads);
  int selectCycleIterative(const int bead0, const int nBisectBeads, Cycle& c);
  void permuteBeads(field<Bead*>& b0, field<Bead*>& b1, Cycle& c);
  void assignParticleLabels();
  void Write();

  std::vector<Bead*> affBeads;
protected:

public:
  // Constructor
  PermBisectIterative(Path &tmpPath, RNG &tmpRNG, vector<Action*> &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void MakeMove();
};


#endif
