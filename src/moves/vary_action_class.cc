#include "vary_action_class.h"

// Accept current move
void VaryAction::Accept()
{
  // Call accept for each action
  action->SetParamSet(param_set_new);
  action->Accept();
}

// Reject current move
void VaryAction::Reject()
{
  // Call reject for each action
  action->SetParamSet(param_set_old);
  action->Reject();

}

// VaryActionion Move
bool VaryAction::Attempt()
{
  // Insert dummy particle
  std::vector<std::pair<uint32_t,uint32_t>> particles;
  particles.push_back(std::make_pair(species_i,0));

  // Get old and set new parameter sets
  param_set_old = action->GetParamSet();
  action->SetParamSet(rng.UnifRand(action->GetNumParamSets())-1);
  param_set_new = action->GetParamSet();

  // Calculate action change
  double old_action = 0.;
  double new_action = 0.;
  if (param_set_old != param_set_new) {
    // Old action
    path.SetMode(OLD_MODE);
    action->SetParamSet(param_set_old);
    old_action += action->GetAction(0, path.n_bead-1, particles, 0);

    // New action
    path.SetMode(NEW_MODE);
    action->SetParamSet(param_set_new);
    new_action += action->GetAction(0, path.n_bead-1, particles, 0);
  }

  double current_action_change = new_action - old_action;
  double log_accept_probablity = -current_action_change;

  // Metropolis reject step
  if (log_accept_probablity < log(rng.UnifRand()))
    return 0;

  return 1;
}
