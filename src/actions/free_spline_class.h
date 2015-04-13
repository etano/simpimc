#ifndef SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_
#define SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_

#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>
#include "../config.h"

/// Bare free particle spline class
class FreeSpline
{
private:
  double i_4_lambda_tau; ///< 1/(4\lambda\tau) for splined action
  uint32_t n_d; ///< Number of dimensions for splined action
  UBspline_1d_d* rho_free_r_spline; ///< Holds the splined action

  /// Creates splined action
  void Init(const double L, const uint32_t n_images);
public:
  // Constructor calls Init
  FreeSpline(const double L, const uint32_t n_images, const uint32_t t_n_d, const double t_i_4_lambda_tau)
    : n_d(t_n_d), i_4_lambda_tau(t_i_4_lambda_tau)
  {
    Init(L,n_images);
  }

  /// Returns the value of g_ij
  double GetGij(const vec<double> &r);

  /// Returns the spatial derivative of g_ij
  double GetGijDGijDr(const vec<double> &r, vec<double> &dgij_dr);
};

#endif // SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_
