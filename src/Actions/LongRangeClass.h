#ifndef LongRangeClass_H
#define LongRangeClass_H

#include "ActionClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>
#include <einspline/multi_nubspline.h>
#include <einspline/nubspline.h>

class LongRange : public Action
{
private:
  int nImages;
  int nOrder;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int offsetA, offsetB;
  int maxLevel;
  int nVal, nTau;

  void ReadFile(string fileName);
protected:

public:
  // Constructor
  LongRange(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Data
  NUgrid* grid;
  Tvector taus;
  arma::field<multi_NUBspline_1d_d*> Ukj, dUkjdBeta;

  // Functions
  virtual void Init(Input &in);
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level);
  virtual void Write();

  // Pair actions
  void CalcUdUVsqz(RealType s, RealType q, RealType z, int level, RealType &U, RealType &dU, RealType &V);
  void CalcUdUVr(Tvector& r, Tvector& rP, int level, RealType &U, RealType &dU, RealType &V);
  void CalcUsqz(RealType s, RealType q, RealType z, int level, RealType &U);
  void CalcUr(Tvector& r, Tvector& rP, int level, RealType &U);
};

#endif
