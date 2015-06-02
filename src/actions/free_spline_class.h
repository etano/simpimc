#ifndef SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_
#define SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_

#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

/// Bare free particle spline class
class FreeSpline
{
private:
  double i_4_lambda_tau; ///< 1/(4\lambda\tau) for splined action
  double i_4_lambda_tau_tau; ///< 1/(4\lambda\tau\tau) for splined action
  UBspline_1d_d* image_action_spline; ///< Holds the splined action
  UBspline_1d_d* d_image_action_d_tau_spline; ///< Holds the tau derivative of the splined action

public:
  // Copy constructor
  FreeSpline(FreeSpline&& fs)
    : image_action_spline(std::move(fs.image_action_spline)), d_image_action_d_tau_spline(std::move(fs.d_image_action_d_tau_spline)), i_4_lambda_tau(std::move(fs.i_4_lambda_tau)), i_4_lambda_tau_tau(std::move(fs.i_4_lambda_tau_tau))
  {}

  FreeSpline& operator=(const FreeSpline& fs) = default;

  /// Constructor
  FreeSpline(const double L, const uint32_t n_images, const double lambda, const double tau, const bool use_tau_derivative)
  {
    // Make constants
    i_4_lambda_tau = 1./(4.*lambda*tau);
    i_4_lambda_tau_tau = 1./(4.*lambda*tau*tau);

    // Setup grid
    Ugrid r_grid;
    r_grid.start = -L/2.;
    r_grid.end = L/2.;
    r_grid.num = 10000;
    double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

    // Make image_action, d_image_action_d_tau, d_image_action_d_r
    vec<double> image_action, d_image_action_d_tau, d_image_action_d_r;
    image_action.zeros(r_grid.num);
    if (use_tau_derivative)
      d_image_action_d_tau.zeros(r_grid.num);
    for (uint32_t i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r_2 = r*r;
      double r_2_i_4_lambda_tau = r_2*i_4_lambda_tau;
      for (uint32_t image=1; image<=n_images; image++) {
        double r_p(r+image*L);
        double r_p_2_i_4_lambda_tau = r_p*r_p*i_4_lambda_tau;
        double d_r_2_r_p_2_i_4_lambda_tau = r_2_i_4_lambda_tau - r_p_2_i_4_lambda_tau;
        double exp_part_p = exp(d_r_2_r_p_2_i_4_lambda_tau);
        double r_m(r-image*L);
        double r_m_2_i_4_lambda_tau = r_m*r_m*i_4_lambda_tau;
        double d_r_2_r_m_2_i_4_lambda_tau = r_2_i_4_lambda_tau - r_m_2_i_4_lambda_tau;
        double exp_part_m = exp(d_r_2_r_m_2_i_4_lambda_tau);
        image_action(i) += exp_part_p + exp_part_m;
        if (use_tau_derivative)
          d_image_action_d_tau(i) += (d_r_2_r_p_2_i_4_lambda_tau*exp_part_p + d_r_2_r_m_2_i_4_lambda_tau*exp_part_m)/tau;
      }
      if (use_tau_derivative)
        d_image_action_d_tau(i) = d_image_action_d_tau(i)/(1.+image_action(i));
      image_action(i) = -log1p(image_action(i));
    }
    BCtype_d xBC = {NATURAL, NATURAL};
    image_action_spline = create_UBspline_1d_d(r_grid, xBC, image_action.memptr());
    if (use_tau_derivative)
      d_image_action_d_tau_spline = create_UBspline_1d_d(r_grid, xBC, d_image_action_d_tau.memptr());
  }

  /// Returns rho_free
  double GetRhoFree(const vec<double> &r)
  {
    return exp(GetLogRhoFree(r));
  }

  /// Returns log(rho_free)
  double GetLogRhoFree(const vec<double> &r)
  {
    double tot = 0.;
    for (const auto &r_d : r) {
      double image_action;
      eval_UBspline_1d_d(image_action_spline,r_d,&image_action);
      tot -= image_action;
    }
    return tot - dot(r,r)*i_4_lambda_tau;
  }

  /// Returns d(rho_free)/dtau
  double GetDRhoFreeDTau(const vec<double> &r)
  {
    return GetDLogRhoFreeDTau(r)*GetRhoFree(r);
  }

  /// Returns dlog(rho_free)/dtau
  double GetDLogRhoFreeDTau(const vec<double> &r)
  {
    double tot = 0.;
    for (const auto &r_d : r) {
      double d_image_action_d_tau;
      eval_UBspline_1d_d(d_image_action_d_tau_spline,r_d,&d_image_action_d_tau);
      tot -= d_image_action_d_tau;
    }
    return tot - dot(r,r)*i_4_lambda_tau_tau;
  }

  /// Sets gradient of rho_free at r and returns rho_free
  double GetGradRhoFree(const vec<double> &r, vec<double> &grad_rho_free)
  {
    vec<double> grad_log_rho_free(r.size());
    double log_rho_free = GetGradLogRhoFree(r,grad_log_rho_free);
    double rho_free = exp(log_rho_free);
    grad_rho_free = grad_log_rho_free*rho_free;
    return rho_free;
  }

  /// Returns gradient of log(rho_free) at r
  double GetGradLogRhoFree(const vec<double> &r, vec<double> &grad_log_rho_free)
  {
    double tot = 0.;
    for (uint32_t d_i=0; d_i<r.size(); d_i++) {
      double image_action, d_image_action_d_r;
      eval_UBspline_1d_d_vg(image_action_spline,r(d_i),&image_action,&d_image_action_d_r);
      grad_log_rho_free(d_i) = -d_image_action_d_r;
      tot -= image_action;
    }
    grad_log_rho_free -= r*2.*i_4_lambda_tau;
    return tot - dot(r,r)*i_4_lambda_tau;
  }

};

#endif // SIMPIMC_ACTIONS_FREE_SPLINE_CLASS_H_
