#ifndef SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Pair action class dervied from David Ceperley's squarer
class DavidPairAction : public PairAction
{
private:
  double du_long_k_0; ///< Constant correction to the k components of the beta derivative of the long range pair action
  double du_long_r_0; ///< Constant correction to spatial beta derivative of the long range pair action
  double u_long_k_0; ///< Constant correction to the k components of the long range pair action
  double u_long_r_0; ///< Constant correction to spatial long range pair action
  double v_long_k_0; ///< Constant correction to the k components of the long range pair potential
  double v_long_r_0; ///< Constant correction to spatial long range pair potential
  uint32_t n_tau; ///< Number of values of tau for which the pair action was computed
  uint32_t n_val; ///< Number of values in the pair action expansion
  vec<double> du_long_k; ///< Vector of k components of the beta derivative of the long range pair action
  vec<double> taus; ///< Vector of tau values for which the pair action was computed
  vec<double> u_long_k; ///< Vector of k components of the long range pair action
  vec<double> v_long_k; ///< Vector of k components of the long range pair potential
  field<multi_NUBspline_1d_d*> u_kj; ///< Spline of k,j component in pair action expansion
  field<multi_NUBspline_1d_d*> du_kj_dbeta; ///< Spline of k,j component in beta derivative of pair action expansion
  NUgrid* grid; ///< Grid for the splines

  /// Reads the file storing the pair action
  virtual void ReadFile(std::string file_name);

  /// Calculate the potential
  virtual double CalcV(double r, double r_p, const uint32_t level);

  /// Calculate the long ranged part of the potential
  virtual double CalcVLong();

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the action
  virtual double CalcULong(const uint32_t b0, const uint32_t b1, const uint32_t level);

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the beta derivative of the action
  virtual double CalcdUdBetaLong();

  /// Calculate the gradient of the action
  vec<double> CalcGradientU(const uint32_t b_i, const uint32_t jB, const uint32_t p_i, const uint32_t p_j, const uint32_t level);

  /// Calculate the Laplacian of the action
  double CalcLaplacianU(const uint32_t b_i, const uint32_t jB, const uint32_t p_i, const uint32_t p_j, const uint32_t level);
public:
  /// Constructor calls PairAction class Init
  DavidPairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

};

#endif // SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_
