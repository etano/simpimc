#ifndef SIMPIMC_ACTIONS_KINETIC_CLASS_H_
#define SIMPIMC_ACTIONS_KINETIC_CLASS_H_

#include "action_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class Kinetic : public Action
{
private:
  int n_images;
  uint max_level;
  std::string species;
  uint species_i, n_part;
  double i_4_lambda_tau;

  // Splines
  field<UBspline_1d_d*> rho_free_r_splines;
  UBspline_1d_d* num_sum_r_spline;
  void SetupSpline();

  // Sums over images
  double GetGaussSum(const double r, const double r2i_4_lambda_tau, const uint slice_diff);
  double GetNumSum(const double r, const double r2i_4_lambda_tau);

protected:

public:
  // Constructor
  Kinetic(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual vec<double> GetActionGradient(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual double GetActionLaplacian(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_KINETIC_CLASS_H_
