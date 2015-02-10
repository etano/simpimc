#ifndef KineticClass_H
#define KineticClass_H

#include "ActionClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class Kinetic : public Action
{
private:
  int nImages, maxLevel;
  string species;
  int iSpecies, nPart;
  double i4LambdaTau;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  UBspline_1d_d* num_sum_r_spline;
  void SetupSpline();

  // Sums over images
  double GetGaussSum(const double &r, const int sliceDiff);
  double GetNumSum(const double &r);

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
  virtual double GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual vec<double> GetActionGradient(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual double GetActionLaplacian(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level);
  virtual void Write();
};

#endif
