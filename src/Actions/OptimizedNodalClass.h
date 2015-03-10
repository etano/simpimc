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
  uint iParamSet, iModel;

  // Rho matrix
  virtual double GetGij(const vec<double> &r, const uint sliceDiff);

  // 1/(4\lambda\tau)
  virtual double Geti4LambdaTau(const uint sliceDiff);

  // Splines
  field<UBspline_1d_d*> rho_node_r_splines;
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
  uint GetParamSet() { return iParamSet; };
  uint GetNumParamSets() { return paramSets.size(); };
  void SetParamSet(uint t_iParamSet) { iParamSet = t_iParamSet; };
  void SetRandomParamSet() { SetParamSet(rng.unifRand(paramSets.size())-1); };

};

#endif
