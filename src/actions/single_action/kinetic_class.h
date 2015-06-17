#ifndef SIMPIMC_ACTIONS_KINETIC_CLASS_H_
#define SIMPIMC_ACTIONS_KINETIC_CLASS_H_

#include "single_action_class.h"
#include "../free_spline_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

/// A kinetic action class
class Kinetic : public SingleAction
{
private:
  std::vector<FreeSpline> rho_free_splines; ///< Holds the splined action for every time slice

  /// Creates splined action for all time slices
  void SetupSpline()
  {
    // Create splines
    uint32_t n_spline = path.n_bead/2 + (path.n_bead%2) + 1;
    rho_free_splines.resize(n_spline);
    rho_free_splines[0] = FreeSpline(path.L, n_images, lambda, path.tau, true);
    #pragma omp parallel for
    for (uint32_t spline_i=1; spline_i<n_spline; ++spline_i)
      rho_free_splines[spline_i] = FreeSpline(path.L, n_images, lambda, path.tau*(spline_i+1), false);
  }
public:
  /// Constructor calls Init
  Kinetic(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    SetupSpline();
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta()
  {
    double tot = n_part*path.n_bead*path.n_d/(2.*path.tau); // Constant term
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      for (uint32_t p_i=0; p_i<n_part; p_i++) {
        vec<double> dr(path.Dr(path(species_i,p_i,b_i),path.GetNextBead(path(species_i,p_i,b_i),1)));
        tot += rho_free_splines[0].GetDLogRhoFreeDTau(dr);
      }
    }

    return tot;
  }

  /// Returns the virial contribution of the action
  virtual double VirialEnergy(const double virial_window_size)
  {
    // Constant term
    double tot = n_part*path.n_bead*path.n_d/(2.*virial_window_size*path.tau);

    // Permutation/winding term
    for (uint32_t p_i=0; p_i<n_part; p_i++) {

      // First bead
      std::shared_ptr<Bead> b_0(path(species_i,p_i,0));
      std::shared_ptr<Bead> b_1(path.GetNextBead(b_0,1));
      std::shared_ptr<Bead> b_2(path.GetNextBead(b_1,1));

      vec<double> dr_i_1(path.Dr(b_1,b_0));
      vec<double> dr_i_L(dr_i_1);
      for (uint32_t skip=1; skip<virial_window_size; skip++) {
        dr_i_L += path.Dr(b_2, b_1);
        b_1 = b_2;
        b_2 = path.GetNextBead(b_2,1);
      }
      tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.tau);

      // Other beads
      for (uint32_t b_i=1; b_i<path.n_bead; b_i++) {
        // Subtract first link and add new last link
        dr_i_L -= dr_i_1;
        dr_i_L += path.Dr(b_2,b_1);

        // Move beads one forward
        b_0 = path.GetNextBead(b_0,1);
        dr_i_1 = path.Dr(path.GetNextBead(b_0,1),b_0);
        b_1 = b_2;
        b_2 = path.GetNextBead(b_2,1);

        // Add to total
        tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.tau);
      }

        //// Image tau derivative
        //vec<double> dr_i(path.Dr(b_0, path.GetPrevBead(b_0,1)));
        //for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
        //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
        //  tot -= GetDLogRhoFreeDTau(dr_i(d_i),r2_i_4_lambda_tau) + r2_i_4_lambda_tau/path.tau;
        //}

        //// Image r derivative
        //vec<double> image_action_gradient(path.n_d);
        //for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
        //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
        //  image_action_gradient(d_i) = GetDLogRhoFreeDR(dr_i(d_i),r2_i_4_lambda_tau) + dr_i(d_i)*.2*i_4_lambda_tau;
        //  r2_i_4_lambda_tau = dr_i_1(d_i)*dr_i_1(d_i)*i_4_lambda_tau;
        //  image_action_gradient(d_i) += GetDLogRhoFreeDR(dr_i_1(d_i),r2_i_4_lambda_tau) + dr_i_1(d_i)*.2*i_4_lambda_tau;
        //}
        //tot -= dot(image_action_gradient,path.Dr(b_0,path(species_i,p_i,0)))/(2.*path.tau);

    }
    return tot;
  }

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    double tot = 0.;
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      if (s_i == species_i) {
        std::shared_ptr<Bead> beadA(path(s_i,p_i,b0));
        std::shared_ptr<Bead> beadF(path.GetNextBead(beadA,b1-b0));
        while(beadA != beadF) {
          std::shared_ptr<Bead> beadB(path.GetNextBead(beadA,skip));
          vec<double> dr(path.Dr(beadA,beadB));
          tot -= rho_free_splines[skip-1].GetLogRhoFree(dr);
          beadA = beadB;
        }
      }
    }

    return tot;
  }

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    vec<double> tot;
    tot.zeros(path.n_d);
    std::shared_ptr<Bead> beadA, beadB, beadC, beadF;
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      if (s_i == species_i) {
        double gauss_prod, rho_free, dist;
        beadA = path(species_i,p_i,b0);
        beadF = path.GetNextBead(beadA,b1-b0);
        while(beadA != beadF) {
          beadB = path.GetPrevBead(beadA,skip);
          vec<double> dr(path.Dr(beadB,beadA));
          tot -= dr;
          beadC = path.GetNextBead(beadA,skip);
          dr = path.Dr(beadA,beadC);
          tot += dr;
          beadA = beadC;
        }
      }
    }

    return 2.*i_4_lambda_level_tau*tot;
  }

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    double tot = 0.;
    vec<double> dr(path.n_d);
    std::shared_ptr<Bead> beadA, beadF;
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      if (s_i == species_i) {
        double gauss_prod, rho_free, dist;
        beadA = path(species_i,p_i,b0);
        beadF = path.GetNextBead(beadA,b1-b0);
        while(beadA != beadF) {
          tot += path.n_d*4.*i_4_lambda_level_tau;
          beadA = path.GetNextBead(beadA,skip);
        }
      }
    }

    return tot;
  }

  /// Write information about the action
  virtual void Write() {}

};

#endif // SIMPIMC_ACTIONS_KINETIC_CLASS_H_
