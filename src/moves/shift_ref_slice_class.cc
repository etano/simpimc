#include "shift_ref_slice_class.h"

void ShiftRefSlice::Init(Input &in)
{
  species = in.GetAttribute<std::string>("species");
  path.GetSpeciesInfo(species,species_i);

  // Generate action list
  std::vector<std::string> species_list;
  species_list.push_back(species);
  GenerateActionList(species_list);
}

// Accept current move
void ShiftRefSlice::Accept()
{
  // Set new ref bead
  path.ref_bead = ref_bead_1;

  // Call accept for each action
  for (auto& action: action_list) {
    if (action->type == "Nodal") // fixme: will affect other species nodal actions
      action->Accept();
  }
}

// Reject current move
void ShiftRefSlice::Reject()
{
  // Reset old ref bead
  path.ref_bead = ref_bead_0;

  // Call reject for each action
  for (auto& action: action_list) {
    if (action->type == "Nodal") // fixme: will affect other species nodal actions
      action->Reject();
  }
}

// ShiftRefSliceion Move
bool ShiftRefSlice::Attempt()
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
      path.SetMode(0);
      path.ref_bead = ref_bead_0;
      old_action += action->GetAction(0, path.n_bead-1, particles, 0);

      // New action
      path.SetMode(1);
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
