#ifndef KineticClass_H
#define KineticClass_H

#include "ActionClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class Kinetic : public Action
{
private:
  int nImages;
  uint maxLevel;
  string species;
  uint iSpecies, nPart;
  double i4LambdaTau;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  UBspline_1d_d* num_sum_r_spline;
  void SetupSpline();

  // Sums over images
  double GetGaussSum(const double r, const double r2i4LambdaTau, const uint sliceDiff);
  double GetNumSum(const double r, const double r2i4LambdaTau);

protected:

public:
  // Constructor
  Kinetic(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual vec<double> GetActionGradient(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual double GetActionLaplacian(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level);
  virtual void Write();
};

#endif
