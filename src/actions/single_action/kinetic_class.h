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
  double d_action_d_beta_const; ///< Constant term for thermal energy
  std::vector<FreeSpline> rho_free_splines; ///< Holds the splined action for every time slice

  /// Creates splined action for all time slices
  void SetupSpline()
  {
    // Create splines
    uint32_t n_spline = species->GetNBead()/2 + (species->GetNBead()%2) + 1;
    rho_free_splines.resize(n_spline);
    rho_free_splines[0] = FreeSpline(path.GetL(), n_images, species->GetLambda(), path.GetTau(), true);
    #pragma omp parallel for
    for (uint32_t spline_i=1; spline_i<n_spline; ++spline_i)
      rho_free_splines[spline_i] = FreeSpline(path.GetL(), n_images, species->GetLambda(), path.GetTau()*(spline_i+1), false);
  }
public:
  /// Constructor calls Init
  Kinetic(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    d_action_d_beta_const = species->GetNPart()*species->GetNBead()*path.GetND()/(2.*path.GetTau());
    SetupSpline();
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta()
  {
    double tot(d_action_d_beta_const); // Constant term
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t b_i=0; b_i<species->GetNBead(); b_i++) {
      for (uint32_t p_i=0; p_i<species->GetNPart(); p_i++) {
        tot += rho_free_splines[0].GetDLogRhoFreeDTau(path.Dr(species->GetBead(p_i,b_i),species->GetBead(p_i,b_i)->GetNextBead(1)));
      }
    }

    return tot;
  }

  /// Returns the virial contribution of the action
  virtual double VirialEnergy(const double virial_window_size)
  {
    // Constant term
    double tot(d_action_d_beta_const/virial_window_size);

    // Permutation/winding term
    for (uint32_t p_i=0; p_i<species->GetNPart(); p_i++) {

      // First bead
      std::shared_ptr<Bead> b_0(species->GetBead(p_i,0));
      std::shared_ptr<Bead> b_1(b_0->GetNextBead(1));
      std::shared_ptr<Bead> b_2(b_1->GetNextBead(1));

      vec<double> dr_i_1(path.Dr(b_1,b_0));
      vec<double> dr_i_L(dr_i_1);
      for (uint32_t skip=1; skip<virial_window_size; skip++) {
        dr_i_L += path.Dr(b_2, b_1);
        b_1 = b_2;
        b_2 = b_2->GetNextBead(1);
      }
      tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.GetTau());

      // Other beads
      for (uint32_t b_i=1; b_i<species->GetNBead(); b_i++) {
        // Subtract first link and add new last link
        dr_i_L -= dr_i_1;
        dr_i_L += path.Dr(b_2,b_1);

        // Move beads one forward
        b_0 = b_0->GetNextBead(1);
        dr_i_1 = path.Dr(b_0->GetNextBead(1),b_0);
        b_1 = b_2;
        b_2 = b_2->GetNextBead(1);

        // Add to total
        tot -= dot(dr_i_L,dr_i_1)*i_4_lambda_tau/(virial_window_size*path.GetTau());
      }

        //// Image tau derivative
        //vec<double> dr_i(path.Dr(b_0, path.GetPrevBead(b_0,1)));
        //for (uint32_t d_i=0; d_i<path.GetND(); d_i++) {
        //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
        //  tot -= GetDLogRhoFreeDTau(dr_i(d_i),r2_i_4_lambda_tau) + r2_i_4_lambda_tau/path.GetTau();
        //}

        //// Image r derivative
        //vec<double> image_action_gradient(path.GetND());
        //for (uint32_t d_i=0; d_i<path.GetND(); d_i++) {
        //  double r2_i_4_lambda_tau = dr_i(d_i)*dr_i(d_i)*i_4_lambda_tau;
        //  image_action_gradient(d_i) = GetDLogRhoFreeDR(dr_i(d_i),r2_i_4_lambda_tau) + dr_i(d_i)*.2*i_4_lambda_tau;
        //  r2_i_4_lambda_tau = dr_i_1(d_i)*dr_i_1(d_i)*i_4_lambda_tau;
        //  image_action_gradient(d_i) += GetDLogRhoFreeDR(dr_i_1(d_i),r2_i_4_lambda_tau) + dr_i_1(d_i)*.2*i_4_lambda_tau;
        //}
        //tot -= dot(image_action_gradient,path.Dr(b_0,path(species_i,p_i,0)))/(2.*path.GetTau());

    }
    return tot;
  }

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    double tot = 0.;
    for (auto& p: particles) {
      if (p.first == species) {
        std::shared_ptr<Bead> bead_a(species->GetBead(p.second,b0));
        std::shared_ptr<Bead> bead_b(bead_a->GetNextBead(b1-b0));
        while(bead_a != bead_b) {
          std::shared_ptr<Bead> next_bead_a(bead_a->GetNextBead(skip));
          tot -= rho_free_splines[skip-1].GetLogRhoFree(path.Dr(bead_a,next_bead_a));
          bead_a = next_bead_a;
        }
      }
    }

    return tot;
  }

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    vec<double> tot(zeros<vec<double>>(path.GetND()));
    for (auto& p: particles) {
      if (p.first == species) {
        double gauss_prod, rho_free, dist;
        std::shared_ptr<Bead> bead_a(species->GetBead(p.second,b0));
        std::shared_ptr<Bead> bead_b(bead_a->GetNextBead(b1-b0));
        while(bead_a != bead_b) {
          tot -= path.Dr(bead_a->GetPrevBead(skip),bead_a);
          std::shared_ptr<Bead> next_bead_a = bead_a->GetNextBead(skip);
          tot += path.Dr(bead_a,next_bead_a);
          bead_a = next_bead_a;
        }
      }
    }

    return 2.*i_4_lambda_level_tau*tot;
  }

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;
    double i_4_lambda_level_tau = i_4_lambda_tau/skip;
    double tot = 0.;
    vec<double> dr(path.GetND());
    for (auto& p: particles) {
      if (p.first == species) {
        double gauss_prod, rho_free, dist;
        std::shared_ptr<Bead> bead_a(species->GetBead(p.second,b0));
        std::shared_ptr<Bead> bead_b(bead_a->GetNextBead(b1-b0));
        while(bead_a != bead_b) {
          tot += path.GetND()*4.*i_4_lambda_level_tau;
          bead_a = bead_a->GetNextBead(skip);
        }
      }
    }

    return tot;
  }

  /// Write information about the action
  virtual void Write() {}

};

#endif // SIMPIMC_ACTIONS_KINETIC_CLASS_H_
