#include "BisectClass.h"

void Bisect::Write()
{
}

void Bisect::MakeMove()
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    nAccept += DoBisect(iPart);
    nAttempt++;
  }
}

// Bisection Move
int Bisect::DoBisect( const int iPart )
{
  unsigned int bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  unsigned int nLevel = int(stepSize);
  //nLevel = 1;
  unsigned int nBisectBeads = pow(2,nLevel); // Number of beads in bisection
  unsigned int bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle
  vector<int> particles;
  particles.push_back(iPart);

  // Set up pointers
  Bead *beadI = path(iPart,bead0);
  Bead *beadF = beadI -> nextB(nBisectBeads);

  // Store coordinates
  Bead *beadA;
  affBeads.clear();
  for(beadA = beadI; beadA != beadF; beadA = beadA -> next) {
    affBeads.push_back(beadA);
    beadA->storeR();
  }

  // Perform the bisection (move exactly through kinetic action)
  int skip;
  double tauEff, sigma, sigma2;
  double lambda = beadI->species.lambda;
  double VA[nLevel], VB[nLevel], dA[nLevel+1], dAold;
  dA[nLevel] = 0.0;
  dAold = 0.0;
  Bead *beadB, *beadC;
  for (int iLevel = nLevel-1; iLevel >= 0; iLevel -= 1) {
    skip = pow(2,iLevel);
    tauEff = path.tau*skip;
    sigma2 = lambda*tauEff;
    sigma = sqrt(sigma2);
    VA[iLevel] = 0.0;
    VB[iLevel] = 0.0;

    beadA = beadI;
    while(beadA != beadF) {
      beadB = beadA -> nextB(skip);
      beadC = beadB -> nextB(skip);

      for (int i=0; i<actionList.size(); ++i)
        VA[iLevel] += actionList[i]->GetAction(beadB->b,beadB->b+skip,particles,iLevel);

      Tvector ac = beadC -> r - beadA -> r;
      path.PutInBox(ac);
      Tvector dr(path.nD);
      rng.normRand(dr, 0, sigma);
      path.PutInBox(dr);
      beadB -> r = beadA -> r + 0.5*ac + dr;

      for (int i=0; i<actionList.size(); ++i)
        VB[iLevel] += actionList[i]->GetAction(beadB->b,beadB->b+skip,particles,iLevel);

      beadA = beadC;
    }

    dA[iLevel] = VB[iLevel] - VA[iLevel];
    dAold = 0.5*(dAold + dA[iLevel+1]);
    // Decide whether or not to accept move (note: only need potential action)
    if ((-dA[iLevel] + dAold) < log(rng.unifRand()))  {
      // Restore coordinates
      path.restoreR(affBeads);
      return 0;
    }
  }

  // Move Accepted
  return 1;
}
