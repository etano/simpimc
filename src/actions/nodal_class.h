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
  std::vector<bool> is_first_time;

  // Nodal distance
  bool use_nodal_distance;
  int dist_type, n_dist_steps;
  double dist_tolerance;
  vec<double> dist, dist_c;
  double GetNodalDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);
  double HybridDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);
  double LineSearchDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);
  double MaxDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);
  double NewtonRaphsonDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  // Nodal action
  vec<double> rho_f, rho_f_c;
  field<vec<double>> grad_rho_f, grad_rho_f_c;
  void SetRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);
  void SetRhoFGradRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  // G_ij matrix
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff) = 0;
  virtual double GetGijDGijDr(const vec<double> &r, const uint32_t slice_diff, vec<double> &dgij_dr) = 0;

  // Actions
  double SimpleAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all);
  double DistanceAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all);
  double DDistanceActionDBeta(const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b);

  // Splines
  virtual void SetupSpline() = 0;

  // RNG
  RNG &rng;

  // Compute magnitude from field of vectors
  template<typename T>
  T FieldVecMag(const field<vec<T>> &field_vec)
  {
    T tot = 0.;
    for (uint32_t i=0; i<field_vec.size(); ++i)
      tot += dot(field_vec(i),field_vec(i));
    return sqrt(tot);
  }

  // Compute magnitude from field of vectors
  template<typename T>
  T FieldVecMag(const field<vec<T>> &field_vec, const uint32_t row)
  {
    T tot = 0.;
    for (uint32_t i=0; i<field_vec.n_cols; ++i)
      tot += dot(field_vec(row,i),field_vec(row,i));
    return sqrt(tot);
  }

public:
  // Constructor
  Nodal(Path &path, RNG &t_rng, Input &in, IO &out)
    : Action(path,in,out), rng(t_rng)
  {}

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);
  virtual double ImportanceWeight();
  virtual void Write() {};
  virtual void Accept();
  virtual void Reject();

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual uint32_t GetParamSet() {};
  virtual uint32_t GetNumParamSets() {};
  virtual void SetParamSet(uint32_t t_param_set_i) {};
  virtual void SetRandomParamSet() {};

};

#endif // SIMPIMC_ACTIONS_NODAL_CLASS_H_
