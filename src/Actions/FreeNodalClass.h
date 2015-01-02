#ifndef FreeNodalClass_H
#define FreeNodalClass_H

#include "NodalClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class FreeNodal : public Nodal
{
private:
  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  void SetupSpline();

protected:
  // Rho matrix
  virtual double GetGij(vec<double> &r, int sliceDiff);

public:
  // Constructor
  FreeNodal(Path &path, RNG &rng, Input &in, IOClass &out)
    : Nodal(path,rng,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);

};

#endif
