#ifndef BarePairActionClass_H
#define BarePairActionClass_H

#include "PairActionClass.h"

class BarePairAction : public PairAction
{
private:

public:
  // Constructor
  BarePairAction(Path &path, Input &in, IOClass &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Parameters
  RealType kCutoff;

  // Splines
  NUBspline_1d_d *v_r_spline, *vLong_r_spline;

  // K values
  Tvector vLong_k;

  // Constant corrections
  RealType vLong_r0, vLong_k0;

  // Grid limits
  RealType r_v_min, r_v_max, r_vLong_min, r_vLong_max;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual RealType CalcV(RealType &r, RealType &rP, int level);
  virtual RealType CalcVLong();
  virtual RealType CalcU(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcULong(int b0, int b1, int level);
  virtual RealType CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcdUdBetaLong();

};

#endif
