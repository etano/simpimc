#ifndef SIMPIMC_ACTIONS_KINETIC_CLASS_H_
#define SIMPIMC_ACTIONS_KINETIC_CLASS_H_

#include "single_action_class.h"
#include "free_spline_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

/// A kinetic action class
class Kinetic : public SingleAction
{
private:
  std::vector<FreeSpline*> rho_free_splines; ///< Holds the splined action for every time slice

  /// Creates splined action for all time slices
  void SetupSpline();
public:
  /// Constructor calls Init
  Kinetic(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    std::cout << "Setting up kinetic action for " << species << "..." << std::endl;
    SetupSpline();

    vec<double> dr(3);
    dr(0) = 1.0;
    dr(1) = 1.0;
    dr(2) = 1.0;
    std::cout << rho_free_splines[0]->GetLogRhoFree(dr);
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the virial contribution of the action
  virtual double VirialEnergy(const double virial_window_size);

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Write information about the action
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_KINETIC_CLASS_H_
