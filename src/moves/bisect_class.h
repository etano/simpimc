#ifndef SIMPIMC_MOVES_BISECT_CLASS_H_
#define SIMPIMC_MOVES_BISECT_CLASS_H_

#include "move_class.h"

class Bisect : public Move
{
private:
  std::string species;
  bool roll_over, adaptive;
  double target_ratio;
  int n_images;
  uint32_t species_i;
  uint32_t n_level, n_bisect_beads;
  uint32_t bead0, bead1;
  double i4_lambda_tau_n_bisect_beads, lambda;
  int ref_accept, ref_attempt;

  std::vector<std::shared_ptr<Bead>> affected_beads;
protected:

public:
  // Constructor
  Bisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();
  virtual void Reset();
  virtual void Write();

};

#endif // SIMPIMC_MOVES_BISECT_CLASS_H_
