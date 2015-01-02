#ifndef PermBisectIterativeClass_H
#define PermBisectIterativeClass_H

#include "MoveClass.h"

class PermBisectIterative : public Move
{
private:
  string species;
  int iSpecies;
  int nImages;
  unsigned int nLevel, nBisectBeads, nPart, nPermPart, nPermType;
  unsigned int bead0, bead1;
  double lambda, i4LambdaTauNBisectBeads, epsilon, logEpsilon, targetRatio;
  bool adaptive, rollOver, fixedNode;
  int refAccept, refAttempt;

  struct Cycle
  {
    double weight, contribution;
    int index, type;
    vec<int> perm, iPerm, part;
  };
  vector<Cycle*> cycles;
  field<Cycle> all_cycles;
  mat<double> t;

  vec<int> permAttempt, permAccept;

  void updatePermTable();
  int selectCycleIterative(Cycle& c);
  void permuteBeads(field< std::shared_ptr<Bead> > &b0, field< std::shared_ptr<Bead> > &b1, Cycle& c);
  void assignParticleLabels();
  void Write();

  std::vector< std::shared_ptr<Bead> > affBeads;
protected:

public:
  // Constructor
  PermBisectIterative(Path &tmpPath, RNG &tmpRNG, std::vector< std::shared_ptr<Action> > &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual int Attempt();
  virtual void Accept();
  virtual void Reject();
  virtual void Reset();
};


#endif
