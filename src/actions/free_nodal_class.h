#ifndef SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_

#include "nodal_class.h"
#include "free_spline_class.h"

/// Bare free particle nodal action class
class FreeNodal : public Nodal
{
private:
  std::vector<FreeSpline> rho_free_splines; ///< Holds the splined action for every time slice

  /// Creates splined action for all time slices
  virtual void SetupSpline();

  /// Returns the value of g_ij
  virtual double GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff);

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double> &dgij_dr);
public:
  // Constructor calls Init
  FreeNodal(Path &path, Input &in, IO &out)
    : Nodal(path,in,out)
  {
    Init(in);
  }

};

#endif // SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
