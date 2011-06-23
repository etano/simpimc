#include "RelabelClass.h"

void Relabel::MakeMove()
{
  nAccept += DoRelabel();
  nAttempt++;
}

// Relabel all the beads, shifting the reference slice
int Relabel::DoRelabel()
{
  path.allAction(path.storeRp);
  path.allAction(path.storeNodeDistancep);

  int bead0 = rng.unifRand(path.nBead) - 1;
  Bead *b;
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    b = path.bead(iPart,bead0);
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      path.bead(iPart,iBead) -> r = b -> rC;
      path.bead(iPart,iBead) -> nDist = b -> nDistC;
      b = b -> next;
    }
  }
  
  if(path.fermi) {
    if(!path.checkConstraint(0)) {
      path.allAction(path.restoreRp);
      path.allAction(path.restoreNodeDistancep);
      return 0;
    }
  }

  return 1;
}
