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

protected:
  int nImages, nOrder, nVal, nTau, maxLevel;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB, offsetA, offsetB;
  bool isConstant, isFirstTime;
  RealType dUdB;

  // Long Range
  int useLongRange;
  RealType kCut;

  // Helpers
  virtual void ReadFile(string fileName) = 0;
  inline void GetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rPrime, NUgrid *g);
  inline void GetConstants(Tvector &rVec, Tvector &rPVec, RealType &r, RealType &rP, RealType &q, RealType &s, RealType &z, RealType &x, RealType &y);

  // Pair actions
  virtual RealType CalcV(Tvector& rVec, Tvector& rPVec, int level) = 0;
  virtual RealType CalcVLong() = 0;
  virtual RealType CalcU(RealType &r, RealType &rP, RealType &s, int level) = 0;
  virtual RealType CalcULong(int b0, int b1, vector<int> &particles, int level) = 0;
  virtual RealType CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level) = 0;
  virtual RealType CalcdUdBetaLong() = 0;

public:
  // Constructor
  PairAction(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {}

  // Functions
  virtual void Init(Input &in);
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level);
  virtual void Write();
  virtual void Accept();
  virtual void Reject();

};

inline void PairAction::GetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rP, NUgrid *g)
{
  rMax = g->end;
  if (r > rMax)
    r = rMax;
  if (rP > rMax)
    rP = rMax;
  rMin = g->start;
  if (rP < rMin)
    rP = rMin;
  if(r < rMin)
    r = rMin;
}

inline void PairAction::GetConstants(Tvector &rVec, Tvector &rPVec, RealType &r, RealType &rP, RealType &q, RealType &s, RealType &z, RealType &x, RealType &y)
{
  r = mag(rVec);
  rP = mag(rPVec);
  q = (r + rP)/2.;
  Tvector dr;
  path.Dr(rVec,rPVec,dr);
  s = mag(dr);
  z = r - rP;
  x = q + 0.5*s;
  y = q - 0.5*s;
}


#endif
