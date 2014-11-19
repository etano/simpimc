#include "DisplaceParticleClass.h"

void DisplaceParticle::Init(Input &in)
{
  species = in.getAttribute<string>("species");
  path.GetSpeciesInfo(species,iSpecies);
  stepSize = in.getAttribute<RealType>("stepSize");
}

// Accept current move
void DisplaceParticle::Accept()
{
  // Move Accepted, so copy new coordinates
  path.storeR(affBeads);
  path.storeRhoKP(affBeads);
  for (int iB=0; iB<path.nBead; ++iB)
    path.storeRhoK(iB,iSpecies);

  // Call accept for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Accept();
}

// Reject current move
void DisplaceParticle::Reject()
{
  // Move rejected, so return old coordinates
  path.restoreR(affBeads);
  path.restoreRhoKP(affBeads);
  for (int iB=0; iB<path.nBead; ++iB)
    path.restoreRhoK(iB,iSpecies);

  // Call reject for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Reject();
}

// Displace particle Move
int DisplaceParticle::Attempt()
{
  // Set which particles are affected by the move
  unsigned int iP = rng.unifRand(path.speciesList[iSpecies]->nPart) - 1;  // Pick particle at random
  vector< pair<int,int> > particles;
  pair<int,int> particle(iSpecies,iP);
  particles.push_back(particle);

  // New sampling
  path.SetMode(1);
  Tvector dr(path.nD);
  rng.unifRand(dr, stepSize);

  // Set which beads are affected by the move
  // and move them
  affBeads.clear();
  for(unsigned int iB=0; iB<path.nBead; ++iB) {
    affBeads.push_back(path(iSpecies,iP,iB));
    path(iSpecies,iP,iB)->move(dr);
  }

  // Calculate action change
  RealType oldAction = 0.;
  RealType newAction = 0.;
  for (int iAction=0; iAction<actionList.size(); ++iAction) {
    // Old action
    path.SetMode(0);
    oldAction += actionList[iAction]->GetAction(0, path.nBead, particles, 0);

    // New action
    path.SetMode(1);
    newAction += actionList[iAction]->GetAction(0, path.nBead, particles, 0);
  }

  RealType logAcceptProb = oldAction - newAction;

  // Metropolis reject step
  if (logAcceptProb < log(rng.unifRand()))
    return 0;
  else
    return 1;
}
