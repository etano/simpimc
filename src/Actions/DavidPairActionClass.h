#ifndef DavidPairActionClass_H
#define DavidPairActionClass_H

#include "PairActionClass.h"

class DavidPairAction : public PairAction
{
private:

public:
  // Constructor
  DavidPairAction(Path &path, Input &in, IOClass &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Data
  NUgrid* grid;
  vec<double> taus;
  field<multi_NUBspline_1d_d*> Ukj, dUkjdBeta;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual double CalcV(double r, double rP, const uint level);
  virtual double CalcVLong() { return 0.; };
  virtual double CalcU(double r, double rP, double s, const uint level);
  virtual double CalcULong(const uint b0, const uint b1, const uint level) { return 0.; };
  virtual double CalcdUdBeta(double r, double rP, double s, const uint level);
  virtual double CalcdUdBetaLong() { return 0.; };

};

#endif
