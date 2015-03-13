#ifndef SIMPIMC_ACTIONS_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_NODAL_CLASS_H_

#include "action_class.h"

class Nodal : public Action
{
private:

protected:
  int n_images;
  uint32_t max_level;
  std::string species;
  uint32_t species_i, n_part;
  double i_4_lambda_tau;
  uint32_t start_b, end_b;

  // Rho matrix
  field<double> rho_f, rho_f_c;
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff) = 0;

  // RNG
  RNG &rng;

public:
  // Constructor
  Nodal(Path &path, RNG &t_rng, Input &in, IO &out)
    : Action(path,in,out), rng(t_rng)
  {}

  // Functions
  virtual void Init(Input &in) {};
  virtual double DActionDBeta();
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);
  virtual void Write() {};
  virtual void Accept();

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual uint32_t GetParamSet() {};
  virtual uint32_t GetNumParamSets() {};
  virtual void SetParamSet(uint32_t t_param_set_i) {};
  virtual void SetRandomParamSet() {};

};

#endif // SIMPIMC_ACTIONS_NODAL_CLASS_H_
