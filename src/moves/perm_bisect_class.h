#ifndef SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
#define simpimc_moves_perm_bisect_class_h_

#include "move_class.h"

class PermBisect : public Move
{
private:
  std::string species;
  int n_images;
  uint32_t species_i;
  uint32_t n_level, n_bisect_beads, n_part, n_perm_part, n_perm_type;
  uint32_t bead0, bead1;
  double lambda, i4_lambda_tau_n_bisect_beads, epsilon, log_epsilon;

  struct Cycle
  {
    double weight, contribution;
    uint32_t index, type;
    vec<uint32_t> perm, i_perm, part;
  };
  std::vector<Cycle*> cycles;
  field<Cycle> all_cycles;
  mat<double> t;

  uint32_t perm_type;
  vec<uint32_t> perm_attempt, perm_accept;

  double ConstructPermTable();
  void UpdatePermTable();
  void BuildCycles();
  uint32_t SelectCycle(const double permTot);
  void PermuteBeads(field<std::shared_ptr<Bead>> &b0, field<std::shared_ptr<Bead>> &b1, const Cycle* const c);
  void AssignParticleLabels();
  void Write();

  std::vector<std::shared_ptr<Bead>> affected_beads;
protected:

public:
  // Constructor
  PermBisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();
};


#endif // SIMPIMC_MOVES_MOVE_CLASS_H_
