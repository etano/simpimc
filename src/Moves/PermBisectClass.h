#ifndef PermBisectClass_H
#define PermBisectClass_H

#include "MoveClass.h"

class PermBisect : public Move
{
private:
  string species;
  int nImages;
  uint iSpecies;
  uint nLevel, nBisectBeads, nPart, nPermPart, nPermType;
  uint bead0, bead1;
  double lambda, i4LambdaTauNBisectBeads, epsilon, logEpsilon;

  struct Cycle
  {
    double weight, contribution;
    uint index, type;
    vec<uint> perm, iPerm, part;
  };
  vector<Cycle*> cycles;
  field<Cycle> all_cycles;
  mat<double> t;

  uint permType;
  vec<uint> permAttempt, permAccept;

  double constructPermTable();
  void updatePermTable();
  void BuildCycles();
  uint selectCycle(const double permTot);
  void permuteBeads(field<std::shared_ptr<Bead> >& b0, field<std::shared_ptr<Bead>>& b1, const Cycle* const c);
  void assignParticleLabels();
  void Write();

  std::vector<std::shared_ptr<Bead> > affBeads;
protected:

public:
  // Constructor
  PermBisect(Path &tmpPath, RNG &tmpRNG, std::vector<std::shared_ptr<Action> > &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();
};


#endif
