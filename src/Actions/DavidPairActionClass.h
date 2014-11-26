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
  virtual double CalcV(double &r, double &rP, int level);
  virtual double CalcVLong() { return 0.; };
  virtual double CalcU(double &r, double &rP, double &s, int level);
  virtual double CalcULong(int b0, int b1, int level) { return 0.; };
  virtual double CalcdUdBeta(double &r, double &rP, double &s, int level);
  virtual double CalcdUdBetaLong() { return 0.; };

};

#endif
