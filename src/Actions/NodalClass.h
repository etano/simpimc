#ifndef NodalClass_H
#define NodalClass_H

#include "ActionClass.h"
#include "../Utils/RNG/RNGClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class Nodal : public Action
{
private:

protected:
  int nImages, maxLevel;
  string species;
  int iSpecies, nPart;
  RealType i4LambdaTau;
  int startB, endB;

  // Rho matrix
  field<RealType> rho_F, rho_F_c;
  RealType GetGij(vec<RealType> &r, int sliceDiff);

  // Ref beads
  vector<Bead*> refBeads;
  field<Bead*> otherBeads;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  void SetupSpline();

  // RNG
  RNG &rng;

public:
  // Constructor
  Nodal(Path &path, RNG &t_rng, Input &in, IOClass &out)
    : Action(path,in,out), rng(t_rng)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector< pair<int,int> > &particles, int level);
  virtual void Write();
  virtual void Accept();

};

#endif
