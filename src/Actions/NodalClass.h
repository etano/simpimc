#ifndef NodalClass_H
#define NodalClass_H

#include "ActionClass.h"

class Nodal : public Action
{
private:

protected:
  int nImages;
  uint maxLevel;
  string species;
  uint iSpecies, nPart;
  double i4LambdaTau;
  uint startB, endB;

  // Rho matrix
  field<double> rho_F, rho_F_c;
  virtual double GetGij(const vec<double> &r, const uint sliceDiff) = 0;

  // RNG
  RNG &rng;

public:
  // Constructor
  Nodal(Path &path, RNG &t_rng, Input &in, IOClass &out)
    : Action(path,in,out), rng(t_rng)
  {}

  // Functions
  virtual void Init(Input &in) {};
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual void Write() {};
  virtual void Accept();

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual uint GetParamSet() {};
  virtual uint GetNumParamSets() {};
  virtual void SetParamSet(uint t_iParamSet) {};
  virtual void SetRandomParamSet() {};

};

#endif
