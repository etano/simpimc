#ifndef NodalClass_H
#define NodalClass_H

#include "ActionClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

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
  double GetGij(vec<double> &r, int sliceDiff);

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
  virtual double DActionDBeta();
  virtual double GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual void Write();
  virtual void Accept();

};

#endif
