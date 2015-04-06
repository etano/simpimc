#ifndef SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Pair action class dervied from David Ceperley's squarer
class DavidPairAction : public PairAction
{
private:
  uint32_t n_tau; ///< Number of values of tau for which the pair action was computed
  uint32_t n_val; ///< Number of values in the pair action expansion
  vec<double> taus; ///< Vector of tau values for which the pair action was computed
  field<multi_NUBspline_1d_d*> ukj; ///< Spline of k,j component in pair action expansion
  field<multi_NUBspline_1d_d*> dukj_dbeta; ///< Spline of k,j component in beta derivative of pair action expansion
  NUgrid* grid; ///< Grid for the splines

  /// Reads the file storing the pair action
  virtual void ReadFile(std::string file_name);

  /// Calculate the potential
  virtual double CalcV(double r, double r_p, const uint32_t level);

  /// Calculate the long ranged part of the potential
  // FIXME: Long ranged potential needs to be implemented
  virtual double CalcVLong() { return 0.; };

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the action
  // FIXME: Long ranged action needs to be implemented
  virtual double CalcULong(const uint32_t b0, const uint32_t b1, const uint32_t level) { return 0.; };

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level);

  /// Calculate the long ranged part of the beta derivative of the action
  // FIXME: Long ranged beta derivative needs to be implemented
  virtual double CalcdUdBetaLong() { return 0.; };

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
