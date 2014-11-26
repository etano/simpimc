#include "BisectClass.h"

void Bisect::Init(Input &in)
{
  nLevel = in.getAttribute<int>("nLevel");
  int maxPossibleLevel = floor(log2(path.nBead));
  if (nLevel > maxPossibleLevel)
    cout << "Warning: nLevel > maxPossibleLevel!" << endl;
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;
  species = in.getAttribute<string>("species");
  path.GetSpeciesInfo(species,iSpecies);

  // Compute constants
  lambda = path.speciesList[iSpecies]->lambda;
  nBisectBeads = 1<<nLevel; // Number of beads in bisection
}

// Accept current move
void Bisect::Accept()
{
  // Move Accepted, so copy new coordinates
  path.storeR(affBeads);
  path.storeRhoKP(affBeads);
  for (int iB=bead0+1; iB<bead1; ++iB)
    path.storeRhoK(iB,iSpecies);

  // Call accept for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Accept();
}

// Reject current move
void Bisect::Reject()
{
  // Move rejected, so return old coordinates
  path.restoreR(affBeads);
  path.restoreRhoKP(affBeads);
  for (int iB=bead0+1; iB<bead1; ++iB)
    path.restoreRhoK(iB,iSpecies);

  // Call reject for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Reject();
}

// Bisection Move
int Bisect::Attempt()
{
  unsigned int iP = rng.unifRand(path.speciesList[iSpecies]->nPart) - 1;  // Pick particle at random
  bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  bead1 = bead0 + nBisectBeads;

  // Set up pointers
  Bead *beadI = path(iSpecies,iP,bead0);
  Bead *beadF = beadI;
  affBeads.clear();
  affBeads.push_back(beadI);
  for (int i=0; i<nBisectBeads; ++i) {
    beadF = beadF->next;
    affBeads.push_back(beadF);
  }

  // Set which particles are affected by the move
  vector< pair<int,int> > particles;
  particles.push_back(std::make_pair(iSpecies,iP));
  if (beadF->p != iP)  // fixme: may be overkill
    particles.push_back(std::make_pair(iSpecies,beadF->p));

  // Perform the bisection (move exactly through kinetic action)
  Bead *beadA, *beadB, *beadC;
  double prevActionChange = 0.;
  double prefactorOfSampleProb = 0.;
  vec<double> rBarOld(path.nD), deltaOld(path.nD), rBarNew(path.nD), deltaNew(path.nD);
  for (int iLevel = nLevel-1; iLevel >= 0; iLevel -= 1) {
    int skip = 1<<iLevel;
    double levelTau = path.tau*skip;
    double sigma2 = lambda*levelTau;
    double sigma = sqrt(sigma2);

    double oldLogSampleProb = 0.;
    double newLogSampleProb = 0.;
    beadA = beadI;
    while(beadA != beadF) {
      // Set beads
      beadB = beadA->nextB(skip);
      beadC = beadB->nextB(skip);

      // Old sampling
      path.SetMode(0);
      path.RBar(beadC, beadA, rBarOld);
      path.Dr(beadB, rBarOld, deltaOld);

      // New sampling
      path.SetMode(1);
      path.RBar(beadC, beadA, rBarNew);
      rng.normRand(deltaNew, 0., sigma);
      path.PutInBox(deltaNew);
      beadB->r = rBarNew + deltaNew;

      // Get sampling probs
      double gaussProdOld = 1.;
      double gaussProdNew = 1.;
      for (int iD=0; iD<path.nD; iD++) {
        double gaussSumOld = 0.;
        double gaussSumNew = 0.;
        for (int image=-nImages; image<=nImages; image++) {
          double distOld = deltaOld(iD) + (double)image*path.L;
          double distNew = deltaNew(iD) + (double)image*path.L;
          gaussSumOld += path.fexp(-0.5*distOld*distOld/sigma2);
          gaussSumNew += path.fexp(-0.5*distNew*distNew/sigma2);
        }
        gaussProdOld *= gaussSumOld;
        gaussProdNew *= gaussSumNew;
      }
      oldLogSampleProb += prefactorOfSampleProb + log(gaussProdOld);
      newLogSampleProb += prefactorOfSampleProb + log(gaussProdNew);

      beadA = beadC;
    }

    // Calculate action change
    double oldAction = 0.;
    double newAction = 0.;
    for (int iAction=0; iAction<actionList.size(); ++iAction) {
      // Old action
      path.SetMode(0);
      oldAction += actionList[iAction]->GetAction(bead0, bead1, particles, iLevel);

      // New action
      path.SetMode(1);
      newAction += actionList[iAction]->GetAction(bead0, bead1, particles, iLevel);
    }

    double logSampleRatio = -newLogSampleProb + oldLogSampleProb;
    double currActionChange = newAction - oldAction;
    double logAcceptProb = logSampleRatio - currActionChange + prevActionChange;

    // Metropolis reject step
    if (logAcceptProb < log(rng.unifRand()))
      return 0;

    prevActionChange = currActionChange;
  }

  return 1;
}
