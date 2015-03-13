#include "vary_optimized_nodal_class.h"

void VaryOptimizedNodal::Init(Input &in)
{
  std::string action_name = in.GetAttribute<std::string>("action_name");

  // Select action from list
  for (auto& t_action : full_action_list) {
    if (t_action->name == action_name)
      action_list.push_back(t_action);
  }

  // Set species
  std::string species = action_list[0]->species_list[0]; // FIXME: Not very general
  path.GetSpeciesInfo(species,species_i);
}

// Accept current move
void VaryOptimizedNodal::Accept()
{
  // Call accept for each action
  for (auto& action: action_list) {
    action->SetParamSet(param_set_1);
    action->Accept();
  }
}

// Reject current move
void VaryOptimizedNodal::Reject()
{
  // Call reject for each action
  for (auto& action: action_list) {
    action->SetParamSet(param_set_0);
    action->Reject();
  }
}

// VaryOptimizedNodalion Move
bool VaryOptimizedNodal::Attempt()
{
  // Insert dummy particle
  std::vector<std::pair<uint32_t,uint32_t>> particles;
  particles.push_back(std::make_pair(species_i,0));

  // Calculate action change
  double old_action = 0.;
  double new_action = 0.;
  for (auto& action: action_list) {
    // Old action
    path.SetMode(0);
    param_set_0 = action->GetParamSet();
    old_action += action->GetAction(0, path.n_bead-1, particles, 0);

    // New action
    path.SetMode(1);
    action->SetRandomParamSet();
    param_set_1 = action->GetParamSet();
    new_action += action->GetAction(0, path.n_bead-1, particles, 0);
  }

  double current_action_change = new_action - old_action;
  double log_accept_probablity = -current_action_change;

  // Metropolis reject step
  if (log_accept_probablity < log(rng.UnifRand()))
    return 0;

  return 1;
}
