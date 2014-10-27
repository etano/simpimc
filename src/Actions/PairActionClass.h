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
  RealType dUdBConstant, VConstant;

  // Long Range
  int useLongRange;
  RealType kCut;

  // Helpers
  virtual void ReadFile(string fileName) = 0;
  inline void GetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rPrime, NUgrid *g);
  inline void SetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rPrime);

  // Pair actions
  virtual RealType CalcV(RealType &r, RealType &rP, int level) = 0;
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
  virtual RealType Potential();
  virtual void Write();
  virtual void Accept();
  virtual void Reject();

};

inline void PairAction::GetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rP, NUgrid *g)
{
  rMin = g->start;
  rMax = g->end;
  SetLimits(rMin,rMax,r,rP);
}

inline void PairAction::SetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rP)
{
  if (r > rMax)
    r = rMax;
  else if (r < rMin)
    r = rMin;

  if (rP > rMax)
    rP = rMax;
  else if (rP < rMin)
    rP = rMin;
}

#endif
