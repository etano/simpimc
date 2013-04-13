#include "DisplaceParticleClass.h"

void DisplaceParticle::MakeMove()
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    nAccept += DoDisplaceParticle(iPart);
    nAttempt++;
  }
}

// Make a single particle displace move
int DisplaceParticle::DoDisplaceParticle( const int iPart )
{
  path.partAction(iPart, path.storeRp); // Store coordinates

  path.mode = 0;
  double V0 = path.getV(iPart);  // Calculate old potential action
  double K0 = path.getK(iPart, 0) + path.getK(iPart, path.nBead-1);  // Calculate old kinetic action
  double A0 = V0 + K0;
  
  // Make a whole particle move  
  rng.unifRand(dr, stepSize);
  for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) path.bead(iPart,iBead) -> move(dr);

  path.mode = 1;
  double V1 = path.getV(iPart);  // Calculate new potential action 
  double K1 = path.getK(iPart, 0) + path.getK(iPart, path.nBead-1);  // Calculate new kinetic action 
  double A1 = V1 + K1;

  // Decide whether or not to accept move
  if ((A0 - A1) < log(rng.unifRand())) {    
    path.partAction(iPart, path.restoreRp); // Restore coordinates
    return 0;   
  }
   
  return 1; 
}
