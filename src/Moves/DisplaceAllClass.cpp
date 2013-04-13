#include "DisplaceAllClass.h"

void DisplaceAll::MakeMove()
{
  nAccept += DoDisplaceAll();
  nAttempt++;
}

// Make a single particle displace move
int DisplaceAll::DoDisplaceAll()
{
  path.allAction(path.storeRp); // Store coordinates

  path.mode = 1;
  double V0 = path.getV();  // Calculate old potential action

  // Make a full move
  rng.unifRand(dr, stepSize);
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      path.bead(iPart,iBead) -> move(dr);
    }
  }

  path.mode = 1;
  double V1 = path.getV();  // Calculate new potential action 

  // Decide whether or not to accept move
  if ((V0 - V1) < log(rng.unifRand())) {
    path.allAction(path.restoreRp); // Restore coordinates
    return 0;
  }

  return 1;
}
