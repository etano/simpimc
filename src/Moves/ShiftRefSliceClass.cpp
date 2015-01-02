#include "ShiftRefSliceClass.h"

void ShiftRefSlice::Init(Input &in)
{
  species = in.getAttribute<string>("species");
  path.GetSpeciesInfo(species,iSpecies);

  // Generate action list
  std::vector<std::string> speciesList;
  speciesList.push_back(species);
  GenerateActionList(speciesList);
}

// Accept current move
void ShiftRefSlice::Accept()
{
  // Set new ref bead
  path.refBead = refBead1;

  // Call accept for each action
  for (auto& action: actionList) {
    if (action->type == "Nodal") // fixme: will affect other species nodal actions
      action->Accept();
  }
}

// Reject current move
void ShiftRefSlice::Reject()
{
  // Reset old ref bead
  path.refBead = refBead0;

  // Call reject for each action
  for (auto& action: actionList) {
    if (action->type == "Nodal") // fixme: will affect other species nodal actions
      action->Reject();
  }
}

// ShiftRefSliceion Move
int ShiftRefSlice::Attempt()
{
  refBead0 = path.refBead;
  refBead1 = rng.unifRand(path.nBead) - 1;  // Pick new ref bead at random

  // Insert dummy particle
  vector< pair<int,int> > particles;
  particles.push_back(std::make_pair(iSpecies,0));

  // Calculate action change
  double oldAction = 0.;
  double newAction = 0.;
  for (auto& action: actionList) {
    // Only check nodal action change
    if (action->type == "Nodal") {
      // Old action
      path.SetMode(0);
      path.refBead = refBead0;
      oldAction += action->GetAction(0, path.nBead-1, particles, 0);

      // New action
      path.SetMode(1);
      path.refBead = refBead1;
      newAction += action->GetAction(0, path.nBead-1, particles, 0);
    }
  }

  double currActionChange = newAction - oldAction;
  double logAcceptProb = -currActionChange;

  // Metropolis reject step
  if (logAcceptProb < log(rng.unifRand()))
    return 0;

  return 1;
}