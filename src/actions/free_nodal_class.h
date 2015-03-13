#ifndef SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_

#include "nodal_class.h"
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
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff);

public:
  // Constructor
  FreeNodal(Path &path, RNG &rng, Input &in, IO &out)
    : Nodal(path,rng,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);

};

#endif // SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
