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
  Tvector taus;
  field<multi_NUBspline_1d_d*> Ukj, dUkjdBeta;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual RealType CalcV(RealType &r, RealType &rP, int level);
  virtual RealType CalcVLong() { return 0.; };
  virtual RealType CalcU(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcULong(int b0, int b1, int level) { return 0.; };
  virtual RealType CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level);
  virtual RealType CalcdUdBetaLong() { return 0.; };

};

#endif
