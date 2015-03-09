#ifndef NodalClass_H
#define NodalClass_H

#include "ActionClass.h"

class Nodal : public Action
{
private:

protected:
  int nImages, maxLevel;
  string species;
  int iSpecies, nPart;
  double i4LambdaTau;
  int startB, endB;

  // Rho matrix
  field<double> rho_F, rho_F_c;
  virtual double GetGij(const vec<double> &r, const int sliceDiff) = 0;

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
  virtual double GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual void Write() {};
  virtual void Accept();

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual int GetParamSet() {};
  virtual int GetNumParamSets() {};
  virtual void SetParamSet(int t_iParamSet) {};
  virtual void SetRandomParamSet() {};

};

#endif
