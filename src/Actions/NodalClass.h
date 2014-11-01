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
  int iSpecies, offset;
  RealType i4LambdaTau;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;

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

  void SetupSpline();
  RealType GetRhoij(Tvector &r, int sliceDiff);

};

#endif
