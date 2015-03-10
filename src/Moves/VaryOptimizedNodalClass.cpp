#include "VaryOptimizedNodalClass.h"

void VaryOptimizedNodal::Init(Input &in)
{
  string actionName = in.getAttribute<string>("actionName");

  // Select action from list
  for (auto& t_action : fullActionList) {
    if (t_action->name == actionName)
      actionList.push_back(t_action);
  }

  // Set species
  string species = actionList[0]->speciesList[0]; // FIXME: Not very general
  path.GetSpeciesInfo(species,iSpecies);
}

// Accept current move
void VaryOptimizedNodal::Accept()
{
  // Call accept for each action
  for (auto& action: actionList) {
    action->SetParamSet(paramSet1);
    action->Accept();
  }
}

// Reject current move
void VaryOptimizedNodal::Reject()
{
  // Call reject for each action
  for (auto& action: actionList) {
    action->SetParamSet(paramSet0);
    action->Reject();
  }
}

// VaryOptimizedNodalion Move
bool VaryOptimizedNodal::Attempt()
{
  // Insert dummy particle
  vector< pair<uint,uint> > particles;
  particles.push_back(std::make_pair(iSpecies,0));

  // Calculate action change
  double oldAction = 0.;
  double newAction = 0.;
  for (auto& action: actionList) {
    // Old action
    path.SetMode(0);
    paramSet0 = action->GetParamSet();
    oldAction += action->GetAction(0, path.nBead-1, particles, 0);

    // New action
    path.SetMode(1);
    action->SetRandomParamSet();
    paramSet1 = action->GetParamSet();
    newAction += action->GetAction(0, path.nBead-1, particles, 0);
  }

  double currActionChange = newAction - oldAction;
  double logAcceptProb = -currActionChange;

  // Metropolis reject step
  if (logAcceptProb < log(rng.unifRand()))
    return 0;

  return 1;
}
