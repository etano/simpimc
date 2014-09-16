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
  arma::field<multi_NUBspline_1d_d*> Ukj, dUkjdBeta;

  // Functions
  virtual void ReadFile(string fileName);

  // Pair actions
  virtual RealType CalcV(Tvector& rVec, Tvector& rPVec, int level) { return 0.; };
  virtual RealType CalcVLong() { return 0.; };
  virtual RealType CalcU(Tvector& rVec, Tvector& rPVec, int level);
  virtual RealType CalcULong(int b0, int b1, vector<int> &particles, int level) { return 0.; };
  virtual RealType CalcdUdBeta(Tvector& rVec, Tvector& rPVec, int level);
  virtual RealType CalcdUdBetaLong() { return 0.; };

};

#endif
