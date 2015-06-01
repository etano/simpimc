#include "vary_action_class.h"

// Accept current move
void VaryAction::Accept()
{
  if (switch_param_sets) {
    // Call accept for each action
    action->SetParamSet(param_set_new);
    action->Accept();
  } else {
    // Call accept for each action
    action->SetParamSet(param_set_old);
    action->Reject();
  }

  // Iterate counters
  n_param_accept(param_set_old,param_set_new) += 1;
  n_param_attempt(param_set_old,param_set_new) += 1;

}

// Reject current move
void VaryAction::Reject()
{
  // Call reject for each action
  action->SetParamSet(param_set_old);
  action->Reject();

  // Iterate counters
  n_param_attempt(param_set_old,param_set_new) += 1;

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

void VaryAction::Write()
{
  mat<double> param_accept_ratio(n_param_sets,n_param_sets);
  for (uint32_t i=0; i<n_param_sets; i++) {
    for (uint32_t j=0; j<n_param_sets; j++) {
      double ratio = n_param_accept(i,j)/n_param_attempt(i,j);
      if (!std::isnormal(ratio))
        ratio = 0.;
      param_accept_ratio(i,j) = ratio;
    }
  }

  // Write
  if (first_time) {
    out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "n_param_accept", n_param_accept);
    out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "n_param_attempt", n_param_attempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "y", param_accept_ratio);
  } else {
    out.AppendDataSet("/Moves/"+name+"/ParamSets/", "n_param_attempt", n_param_attempt);
    out.AppendDataSet("/Moves/"+name+"/ParamSets/", "n_param_accept", n_param_accept);
    out.AppendDataSet("/Moves/"+name+"/ParamSets/", "y", param_accept_ratio);
  }

  // Reset
  Reset();

  Move::Write();
}

void VaryAction::Reset()
{
  // Reset counters
  n_param_accept.zeros();
  n_param_attempt.zeros();
}
