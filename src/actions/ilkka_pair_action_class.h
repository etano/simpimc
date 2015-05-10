#ifndef SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Pair action class derived from Ilkka Kylanpaa's squarer
class IlkkaPairAction : public PairAction
{
private:
  double du_long_k_0; ///< Constant correction to the k components of the beta derivative of the long range pair action
  double du_long_r_0; ///< Constant correction to spatial beta derivative of the long range pair action
  double r_u_max; ///< Maximum distance possible in the pair action spline
  double r_u_min; ///< Minimum distance possible in the pair action spline
  double r_v_max; ///< Maximum distance possible in the pair potential spline
  double r_v_min; ///< Minimum distance possible in the pair potential spline
  double r_du_max; ///< Maximum distance possible in the beta derivative of the pair action spline
  double r_du_min; ///< Minimum distance possible in the beta derivative of the pair action spline
  double r_v_long_max; ///< Maximum distance possible in the long range pair potential spline
  double r_v_long_min; ///< Minimum distance possible in the long range pair potential spline
  double u_long_k_0; ///< Constant correction to the k components of the long range pair action
  double u_long_r_0; ///< Constant correction to spatial long range pair action
  double v_long_k_0; ///< Constant correction to the k components of the long range pair potential
  double v_long_r_0; ///< Constant correction to spatial long range pair potential
  vec<double> u_long_k; ///< Vector of k components of the long range pair action
  vec<double> du_long_k; ///< Vector of k components of the beta derivative of the long range pair action
  vec<double> v_long_k; ///< Vector of k components of the long range pair potential
  NUBspline_1d_d *u_long_r_spline; ///< Spline of the long range pair action
  NUBspline_1d_d *du_long_r_spline; ///< Spline of the beta derivative of the long range pair action
  NUBspline_1d_d *v_r_spline; ///< Spline of the pair potential
  NUBspline_1d_d *v_long_r_spline; ///< Spline of the long range pair potential
  NUBspline_2d_d *u_xy_spline; ///< Spline of the pair action
  NUBspline_2d_d *du_xy_spline; ///< Spline of the beta derivative of the pair action

  /// Reads the file storing the pair action
  virtual void ReadFile(std::string file_name);

  /// Calculate the potential
  virtual double CalcV(double r, double r_p, const uint32_t level);

  /// Calculate the long ranged part of the potential
  virtual double CalcVLong();

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the action
  virtual double CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level);

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the beta derivative of the action
  virtual double CalcdUdBetaLong();

  /// Calculate the gradient of the action for the particle pair p_i, p_j in the direction of particle p_i
  virtual vec<double> CalcGradientU(const uint32_t b_i, const uint32_t b_j, const uint32_t p_i, const uint32_t p_j, const uint32_t level);

  /// Calculate the gradient of the long ranged part of the action for all particles
  virtual vec<double> CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level);

  /// Calculate the gradient of the long ranged part of the action in the direction of p_i
  virtual vec<double> CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t p_i, const uint32_t level);

public:
  /// Constructor only calls Init
  IlkkaPairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }
};

#endif // SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
