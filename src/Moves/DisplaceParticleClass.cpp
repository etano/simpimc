#include "DisplaceParticleClass.h"

void DisplaceParticle::Init(Input &in)
{
  species = in.getAttribute<string>("species");
  path.GetSpeciesInfo(species,iSpecies);
  stepSize = in.getAttribute<double>("stepSize",path.L/10.);

  // Generate action list
  std::vector<std::string> speciesList;
  speciesList.push_back(species);
  GenerateActionList(speciesList);
}

// Accept current move
void DisplaceParticle::Accept()
{
  // Move Accepted, so copy new coordinates
  path.storeR(affBeads);
  path.storeRhoKP(affBeads);
  for (uint iB=0; iB<path.nBead; ++iB)
    path.storeRhoK(iB,iSpecies);

  // Call accept for each action
  for (auto& action: actionList)
    action->Accept();
}

// Reject current move
void DisplaceParticle::Reject()
{
  // Move rejected, so return old coordinates
  path.restoreR(affBeads);
  path.restoreRhoKP(affBeads);
  for (uint iB=0; iB<path.nBead; ++iB)
    path.restoreRhoK(iB,iSpecies);

  // Call reject for each action
  for (auto& action: actionList)
    action->Reject();
}

// Displace particle Move
bool DisplaceParticle::Attempt()
{
  // Set which particles are affected by the move
  uint iP = rng.unifRand(path.speciesList[iSpecies]->nPart) - 1;  // Pick particle at random
  vector< pair<uint,uint> > particles;
  pair<uint,uint> particle(iSpecies,iP);
  particles.push_back(particle);

  // New sampling
  path.SetMode(1);
  vec<double> dr(path.nD);
  rng.unifRand(dr, stepSize);

  // Set which beads are affected by the move
  // and move them
  affBeads.clear();
  for(uint iB=0; iB<path.nBead; ++iB) {
    affBeads.push_back(path(iSpecies,iP,iB));
    path(iSpecies,iP,iB)->move(dr);
  }

  // Calculate action change
  double oldAction = 0.;
  double newAction = 0.;
  for (auto& action: actionList) {
    // Old action
    path.SetMode(0);
    oldAction += action->GetAction(0, path.nBead, particles, 0);

    // New action
    path.SetMode(1);
    newAction += action->GetAction(0, path.nBead, particles, 0);
  }

  double logAcceptProb = oldAction - newAction;

  // Metropolis reject step
  if (logAcceptProb < log(rng.unifRand()))
    return 0;
  else
    return 1;
}
