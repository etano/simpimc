#ifndef OptimizedNodalClass_H
#define OptimizedNodalClass_H

#include "NodalClass.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class OptimizedNodal : public Nodal
{
private:

protected:
  // Variational parameter sets
  vector< vector<double> > paramSets;
  int iParamSet, iModel;

  // Rho matrix
  virtual double GetGij(vec<double> &r, int sliceDiff);

  // 1/(4\lambda\tau)
  virtual double Geti4LambdaTau(int sliceDiff);

  // Splines
  field<UBspline_1d_d*> rho_node_r2_splines;
  virtual void SetupSpline();

public:
  // Constructor
  OptimizedNodal(Path &path, RNG &rng, Input &in, IOClass &out)
    : Nodal(path,rng,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual void Write();
  int GetParamSet() { return iParamSet; };
  int GetNumParamSets() { return paramSets.size(); };
  void SetParamSet(int t_iParamSet) { iParamSet = t_iParamSet; };
  void SetRandomParamSet() { SetParamSet(rng.unifRand(paramSets.size())-1); };

};

#endif
