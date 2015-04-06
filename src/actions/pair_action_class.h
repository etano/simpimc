#ifndef SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_

#include "action_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>
#include <einspline/multi_nubspline.h>
#include <einspline/nubspline.h>
#include <algorithm>

/// Action class for actions between two species of particles
class PairAction : public Action
{
protected:
  bool is_constant; ///< Whether or not the action is constant
  bool is_first_time; ///< Whether or not this is the first time evaluating the action
  bool use_long_range; ///< Whether or not this is a long ranged action
  double dUdB_constant; ///< Constant used in calculating the beta derivative
  double k_cut; ///< Cutoff in k space for the action
  double potential_constant; ///< Constant used in the calculating the potential
  uint32_t n_order; ///< Order of action expansion if used
  uint32_t species_a_i; ///< Index of species affected by the action
  uint32_t species_b_i; ///< Index of species affected by the action
  std::string species_a; ///< Name of species affected by the action
  std::string species_b; ///< Name of species affected by the action

  /// Reads the file storing the pair action
  virtual void ReadFile(std::string file_name) = 0;

  /// Sets the limits of the pair action and applies them
  inline void GetLimits(double &r_min, double &r_max, double &r, double &r_p, const NUgrid* const g)
  {
    r_min = g->start;
    r_max = g->end;
    SetLimits(r_min,r_max,r,r_p);
  }

  /// Applies the limits of the pair action
  inline void SetLimits(const double r_min, const double r_max, double &r, double &r_p)
  {
    if (r > r_max)
      r = r_max;
    else if (r < r_min)
      r = r_min;

    if (r_p > r_max)
      r_p = r_max;
    else if (r_p < r_min)
      r_p = r_min;
  }

  /// Generate all possible pairs that the action affects
  void GenerateParticlePairs(const std::vector<std::pair<uint32_t,uint32_t>> &particles, std::vector<uint32_t> &particles_a, std::vector<uint32_t> &particles_b, std::vector<std::pair<uint32_t,uint32_t>> &particle_pairs);

  /// Calculate the potential
  virtual double CalcV(double r, double r_p, const uint32_t level) = 0;

  /// Calculate the long ranged part of the potential
  virtual double CalcVLong() = 0;

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level) = 0;

  /// Calculate the long ranged part of the action
  virtual double CalcULong(const uint32_t b0, const uint32_t b1, const uint32_t level) = 0;

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level) = 0;

  /// Calculate the long ranged part of the beta derivative of the action
  virtual double CalcdUdBetaLong() = 0;

  /// Calculate the gradient of the action
  vec<double> CalcGradientU(const uint32_t b_i, const uint32_t jB, const uint32_t p_i, const uint32_t p_j, const uint32_t level);

  /// Calculate the Laplacian of the action
  double CalcLaplacianU(const uint32_t b_i, const uint32_t jB, const uint32_t p_i, const uint32_t p_j, const uint32_t level);

public:
  /// Constructor only instantiates the parent Action class
  PairAction(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {}

  /// Initialize the action
  virtual void Init(Input &in);

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the potential of the action for the whole path
  virtual double Potential();

  /// Returns the importance weight of the action for the whole path
  virtual double ImportanceWeight();

  /// Write information about the action
  virtual void Write();

  /// Accepts relevant information about the action
  virtual void Accept();

  /// Rejects relevant information about the action
  virtual void Reject();
};

#endif // SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_
