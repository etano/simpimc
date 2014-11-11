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
  int iSpecies, offset, nPart;
  RealType i4LambdaTau;
  int startB, endB;

  // Rho matrix
  field<Tmatrix> g, g_c, c, c_c;
  field<bool> c_good;
  field<RealType> rho_F, rho_F_c;
  void SetCij(int iB);
  RealType GetGij(Tvector &r, int sliceDiff);

  // Ref beads
  vector<Bead*> refBeads;
  field<Bead*> otherBeads;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  void SetupSpline();

public:
  // Constructor
  Nodal(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level);
  virtual void Write();
  virtual void Accept();

};

#endif
