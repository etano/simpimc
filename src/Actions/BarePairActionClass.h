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
  double kCutoff;

  // Splines
  NUBspline_1d_d *v_r_spline, *vLong_r_spline;

  // K values
  vec<double> vLong_k;

  // Constant corrections
  double vLong_r0, vLong_k0;

  // Grid limits
  double r_v_min, r_v_max, r_vLong_min, r_vLong_max;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual double CalcV(double &r, double &rP, const int level);
  virtual double CalcVLong();
  virtual double CalcU(double &r, double &rP, double &s, const int level);
  virtual double CalcULong(const int b0, const int b1, const int level);
  virtual double CalcdUdBeta(double &r, double &rP, double &s, const int level);
  virtual double CalcdUdBetaLong();

};

#endif
