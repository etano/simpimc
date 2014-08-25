#ifndef PairActionClass_H
#define PairActionClass_H

#include "ActionClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>
#include <einspline/multi_nubspline.h>
#include <einspline/nubspline.h>
#include <algorithm>

class PairAction : public Action
{
private:
  int nImages, nOrder, nVal, nTau, maxLevel;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB, offsetA, offsetB;

  // Long Range
  int useLongRange;
  RealType kCut;
  Tvector Ukl, dUkl;
  vector<Tvector> ks;
  vector<RealType> magKs;
  vector< pair<RealType,int> > uniqueKs;
  vector<Ivector> kis;

  void ReadFile(string fileName);
protected:

public:
  // Constructor
  PairAction(Path &path, Input &in, IOClass &out)
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
