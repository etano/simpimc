#ifndef SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

class IlkkaPairAction : public PairAction
{
private:

public:
  // Constructor
  IlkkaPairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Parameters
  double k_cutoff;

  // Splines
  NUBspline_1d_d *u_long_r_spline, *du_long_r_spline, *v_r_spline, *v_long_r_spline;
  NUBspline_2d_d *u_xy_spline, *du_xy_spline;

  // K values
  vec<double> u_long_k, du_long_k, v_long_k;

  // Constant corrections
  double u_long_r_0, u_long_k_0, du_long_r_0, du_long_k_0, v_long_r_0, v_long_k_0;

  // Grid limits
  double r_u_min, r_u_max, r_du_min, r_du_max, r_v_min, r_v_max, r_v_long_min, r_v_long_max;

  // Functions
  virtual void ReadFile(std::string file_name);

  // Pair actions
  virtual double CalcV(double r, double r_p, const uint level);
  virtual double CalcVLong();
  virtual double CalcU(double r, double r_p, double s, const uint level);
  virtual double CalcULong(const uint b0, const uint b1, const uint level);
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint level);
  virtual double CalcdUdBetaLong();

};

#endif // SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
