#ifndef ImportancePairActionClass_H
#define ImportancePairActionClass_H

#include "PairActionClass.h"

class ImportancePairAction : public PairAction
{
private:

public:
  // Constructor
  ImportancePairAction(Path &path, Input &in, IOClass &out)
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
  virtual double CalcV(double r, double rP, const uint level);
  virtual double CalcVLong();
  virtual double CalcU(double r, double rP, double s, const uint level);
  virtual double CalcULong(const uint b0, const uint b1, const uint level);
  virtual double CalcdUdBeta(double r, double rP, double s, const uint level);
  virtual double CalcdUdBetaLong();

};

#endif
