#ifndef SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Simple diagonal pair action
class BarePairAction : public PairAction
{
private:

public:
  // Constructor
  BarePairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Parameters
  double k_cutoff;

  // Splines
  NUBspline_1d_d *v_r_spline, *v_long_r_spline;

  // K values
  vec<double> v_long_k;

  // Constant corrections
  double v_long_r_0, v_long_k_0;

  // Grid limits
  double r_v_min, r_v_max, r_v_long_min, r_v_long_max;

  // Functions
  virtual void ReadFile(std::string file_name);

  // Pair actions
  virtual double CalcV(double r, double r_p, const uint32_t level);
  virtual double CalcVLong();
  virtual double CalcU(double r, double r_p, double s, const uint32_t level);
  virtual double CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level);
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level);
  virtual double CalcdUdBetaLong();

};

#endif // SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_
