#ifndef IlkkaPairActionClass_H
#define IlkkaPairActionClass_H

#include "PairActionClass.h"

class IlkkaPairAction : public PairAction
{
private:

public:
  // Constructor
  IlkkaPairAction(Path &path, Input &in, IOClass &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Parameters
  int nOrder;
  RealType kCutoff;

  // Splines
  NUBspline_1d_d *uLong_r_spline, *duLong_r_spline, *v_r_spline, *vLong_r_spline;
  NUBspline_2d_d *u_xy_spline, *du_xy_spline;

  // K values
  Tvector uLong_k, duLong_k, vLong_k;

  // Constant corrections
  RealType uLong_r0, uLong_k0, duLong_r0, duLong_k0, vLong_r0, vLong_k0;

  // Grid limits
  RealType r_u_min, r_u_max, r_du_min, r_du_max, r_v_min, r_v_max, r_vLong_min, r_vLong_max;

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
