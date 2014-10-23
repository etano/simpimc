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
  string type1, type2;
  int nOrder;

  // Read values
  Tvector r_u, x_u, y_u, r_du, x_du, y_du, r_v, r_vLong;
  Tvector uLong_r, uLong_k, duLong_r, duLong_k, v_r, vLong_r, vLong_k;
  Tmatrix u_xy, du_xy;

  // Grids
  NUgrid *r_u_grid, *x_u_grid, *y_u_grid;
  NUgrid *r_du_grid, *x_du_grid, *y_du_grid;
  NUgrid *r_v_grid, *r_vLong_grid;

  // Splines
  NUBspline_1d_d *uLong_r_spline, *duLong_r_spline, *v_r_spline, *vLong_r_spline;
  NUBspline_2d_d *u_xy_spline, *du_xy_spline;

  // Constants
  RealType uLong_r0, uLong_k0, duLong_r0, duLong_k0, vLong_r0, vLong_k0;
  RealType kCutoff;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual RealType CalcV(RealType &r, RealType &rP, int level);
  virtual RealType CalcVLong();
  virtual RealType CalcU(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcULong(int b0, int b1, vector<int> &particles, int level);
  virtual RealType CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcdUdBetaLong();

};

#endif
