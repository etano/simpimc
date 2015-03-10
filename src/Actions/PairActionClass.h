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
  int nImages;
  uint nOrder, nVal, nTau, maxLevel;
  string speciesA, speciesB;
  uint iSpeciesA, iSpeciesB;
  bool isConstant, isFirstTime;
  double dUdBConstant, VConstant;

  // Long Range
  bool useLongRange;
  double kCut;

  // Helpers
  virtual void ReadFile(string fileName) = 0;
  inline void GetLimits(double &rMin, double &rMax, double &r, double &rPrime, const NUgrid* const g);
  inline void SetLimits(const double rMin, const double rMax, double &r, double &rPrime);
  void GenerateParticlePairs(const vector<pair<uint,uint> >& particles, vector<uint> &particlesA, vector<uint> &particlesB, vector< pair<uint,uint> > &particlePairs);

  // Pair actions
  virtual double CalcV(double r, double rP, const uint level) = 0;
  virtual double CalcVLong() = 0;
  virtual double CalcU(double r, double rP, double s, const uint level) = 0;
  virtual double CalcULong(const uint b0, const uint b1, const uint level) = 0;
  virtual double CalcdUdBeta(double r, double rP, double s, const uint level) = 0;
  virtual double CalcdUdBetaLong() = 0;
  vec<double> CalcGradientU(const uint iB, const uint jB, const uint iP, const uint jP, const uint level);
  double CalcLaplacianU(const uint iB, const uint jB, const uint iP, const uint jP, const uint level);

public:
  // Constructor
  PairAction(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {}

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual vec<double> GetActionGradient(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual double GetActionLaplacian(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual double Potential();
  virtual double ImportanceWeight();
  virtual void Write();
  virtual void Accept();
  virtual void Reject();

};

inline void PairAction::GetLimits(double &rMin, double &rMax, double &r, double &rP, const NUgrid* const g)
{
  rMin = g->start;
  rMax = g->end;
  SetLimits(rMin,rMax,r,rP);
}

inline void PairAction::SetLimits(const double rMin, const double rMax, double &r, double &rP)
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
