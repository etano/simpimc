#ifndef SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_

#include "action_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>
#include <einspline/multi_nubspline.h>
#include <einspline/nubspline.h>
#include <algorithm>

class PairAction : public Action
{
private:

protected:
  int n_images;
  uint n_order, n_val, n_tau, max_level;
  std::string species_a, species_b;
  uint species_a_i, species_b_i;
  bool is_constant, is_first_time;
  double dUdB_constant, potential_constant;

  // Long Range
  bool use_long_range;
  double k_cut;

  // Helpers
  virtual void ReadFile(std::string file_name) = 0;
  inline void GetLimits(double &r_min, double &r_max, double &r, double &r_prime, const NUgrid* const g);
  inline void SetLimits(const double r_min, const double r_max, double &r, double &r_prime);
  void GenerateParticlePairs(const std::vector<std::pair<uint,uint>> &particles, std::vector<uint> &particles_a, std::vector<uint> &particles_b, std::vector<std::pair<uint,uint>> &particle_pairs);

  // Pair actions
  virtual double CalcV(double r, double r_p, const uint level) = 0;
  virtual double CalcVLong() = 0;
  virtual double CalcU(double r, double r_p, double s, const uint level) = 0;
  virtual double CalcULong(const uint b0, const uint b1, const uint level) = 0;
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint level) = 0;
  virtual double CalcdUdBetaLong() = 0;
  vec<double> CalcGradientU(const uint b_i, const uint jB, const uint p_i, const uint p_j, const uint level);
  double CalcLaplacianU(const uint b_i, const uint jB, const uint p_i, const uint p_j, const uint level);

public:
  // Constructor
  PairAction(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {}

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual vec<double> GetActionGradient(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual double GetActionLaplacian(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  virtual double Potential();
  virtual double ImportanceWeight();
  virtual void Write();
  virtual void Accept();
  virtual void Reject();

};

inline void PairAction::GetLimits(double &r_min, double &r_max, double &r, double &r_p, const NUgrid* const g)
{
  r_min = g->start;
  r_max = g->end;
  SetLimits(r_min,r_max,r,r_p);
}

inline void PairAction::SetLimits(const double r_min, const double r_max, double &r, double &r_p)
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

#endif // SIMPIMC_ACTIONS_PAIR_ACTION_CLASS_H_
