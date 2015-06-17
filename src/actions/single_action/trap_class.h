#ifndef SIMPIMC_ACTIONS_TRAP_CLASS_H_
#define SIMPIMC_ACTIONS_TRAP_CLASS_H_

#include "single_action_class.h"

/// Harmonic well action class
class Trap : public SingleAction
{
private:
  double omega; ///< Harmonic frequency of the trap
  double cofactor_a; ///< Constant from 4th order expansion of the action
  double cofactor_b; ///< Constant from 4th order expansion of the beta derivative of the action
public:
  /// Constructor calls Init
  Trap(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    omega = in.GetAttribute<double>("omega");
    out.Write("/Actions/"+name+"/omega", omega);

    cofactor_a = 0.5*path.tau*omega*omega*(1. + path.tau*path.tau*omega*omega/12.);
    cofactor_b = 0.5*omega*omega*(1. + 3.*path.tau*path.tau*omega*omega/12.);
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta()
  {
    double tot = 0.;
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t p_i=0; p_i<path.n_part; p_i+=1) {
      for (uint32_t b_i=0; b_i<path.n_bead; b_i+=1) {
        tot += dot(path(species_i,p_i,b_i)->r, path(species_i,p_i,b_i)->r);
      }
    }

    return cofactor_b*tot;
  }

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    if (level > max_level)
      return 0.;

    bool check = false;
    for (uint32_t p=0; p<particles.size(); ++p) {
      if (particles[p].first == species_i) {
        check = true;
        break;
      }
    }
    if (!check)
      return 0.;

    uint32_t skip = 1<<level;
    double tot = 0.;
    for (uint32_t p=0; p<particles.size(); ++p) {
      uint32_t s_i = particles[p].first;
      uint32_t p_i = particles[p].second;
      for (uint32_t b_i=b0; b_i<b1; b_i+=skip) {
        vec<double> dr(path.GetR(path(s_i,p_i,b_i)));
        tot += dot(dr, dr);
      }
    }

    return cofactor_a*skip*tot;
  }

  /// Write information about the action
  virtual void Write() {}

};

#endif // SIMPIMC_ACTIONS_TRAP_CLASS_H_
