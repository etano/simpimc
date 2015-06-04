#ifndef SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
#define SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_

#include "single_species_move_class.h"

/// This move attempts to shift the reference slice
class ShiftRefSlice : public SingleSpeciesMove
{
private:
  uint32_t ref_bead_0; ///< Original refence slice
  uint32_t ref_bead_1; ///< New reference slice

  /// Accept the move
  virtual void Accept()
  {
    // Set new ref bead
    path.ref_bead = ref_bead_1;

    // Call accept for each action
    for (auto& action: action_list) {
      if (action->type == "Nodal") // fixme: will affect other species nodal actions
        action->Accept();
    }
  }

  /// Attempt the move
  virtual bool Attempt()
  {
    ref_bead_0 = path.ref_bead;
    ref_bead_1 = rng.UnifRand(path.n_bead) - 1;  // Pick new ref bead at random

    // Insert dummy particle
    std::vector<std::pair<uint32_t,uint32_t>> particles;
    particles.push_back(std::make_pair(species_i,0));

    // Calculate action change
    double old_action = 0.;
    double new_action = 0.;
    for (auto& action: action_list) {
      // Only check nodal action change
      if (action->type == "Nodal") {
        // Old action
        path.SetMode(OLD_MODE);
        path.ref_bead = ref_bead_0;
        old_action += action->GetAction(0, path.n_bead-1, particles, 0);

        // New action
        path.SetMode(NEW_MODE);
        path.ref_bead = ref_bead_1;
        new_action += action->GetAction(0, path.n_bead-1, particles, 0);
      }
    }

    double current_action_change = new_action - old_action;
    double log_accept_probablity = -current_action_change;

    // Metropolis reject step
    if (log_accept_probablity < log(rng.UnifRand()))
      return 0;

    return 1;
  }

  /// Rejects the move
  virtual void Reject()
  {
    // Reset old ref bead
    path.ref_bead = ref_bead_0;

    // Call reject for each action
    for (auto& action: action_list) {
      if (action->type == "Nodal") // fixme: will affect other species nodal actions
        action->Reject();
    }
  }

public:
  /// Constructor only instantiates parent
  ShiftRefSlice(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {}
};

#endif // SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
