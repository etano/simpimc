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
  int iSpeciesA, iSpeciesB;
  bool isConstant, isFirstTime;
  double dUdBConstant, VConstant;

  // Long Range
  int useLongRange;
  double kCut;

  // Helpers
  virtual void ReadFile(string fileName) = 0;
  inline void GetLimits(double &rMin, double &rMax, double &r, double &rPrime, NUgrid *g);
  inline void SetLimits(double &rMin, double &rMax, double &r, double &rPrime);
  void GenerateParticlePairs(const vector<pair<int,int> > &particles, vector<int> &particlesA, vector<int> &particlesB, vector< pair<int,int> > &particlePairs);

  // Pair actions
  virtual double CalcV(double &r, double &rP, const int level) = 0;
  virtual double CalcVLong() = 0;
  virtual double CalcU(double &r, double &rP, double &s, const int level) = 0;
  virtual double CalcULong(const int b0, const int b1, const int level) = 0;
  virtual double CalcdUdBeta(double &r, double &rP, double &s, const int level) = 0;
  virtual double CalcdUdBetaLong() = 0;
  vec<double> CalcGradientU(const int &iB, const int &jB, const int &iP, const int &jP, const int &level);
  double CalcLaplacianU(const int &iB, const int &jB, const int &iP, const int &jP, const int &level);

public:
  // Constructor
  PairAction(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {}

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual vec<double> GetActionGradient(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual double GetActionLaplacian(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual double Potential();
  virtual double ImportanceWeight();
  virtual void Write();
  virtual void Accept();
  virtual void Reject();

};

inline void PairAction::GetLimits(double &rMin, double &rMax, double &r, double &rP, NUgrid *g)
{
  rMin = g->start;
  rMax = g->end;
  SetLimits(rMin,rMax,r,rP);
}

inline void PairAction::SetLimits(double &rMin, double &rMax, double &r, double &rP)
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
