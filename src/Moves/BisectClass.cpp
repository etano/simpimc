#include "BisectClass.h"

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
  unsigned int nBisectBeads = pow(2,nLevel); // Number of beads in bisection
  unsigned int bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle  
  
  // Set up pointers
  Bead *beadI = path.bead(iPart,bead0);
  Bead *beadF = beadI -> nextB(nBisectBeads);

  // Store coordinates   
  Bead *beadA;
  affBeads.clear();
  for(beadA = beadI; beadA != beadF; beadA = beadA -> next) {
    affBeads.push_back(beadA);
    beadA -> storeR();   
  }

  // Figure out which time slices' nodes are affected
  unsigned int nodeBead0, nodeBead1;
  if(rollOver) {
    nodeBead0 = 1;
    nodeBead1 = path.nBead;
  } else {
    nodeBead0 = bead0 + 1;
    nodeBead1 = bead1;
  }
  
  // Old nodal action
  double N0(0.0);
  if(path.useNodeDist) {
    for (unsigned int iBead = nodeBead0; iBead < nodeBead1; iBead += 1) {
      for (unsigned int jPart = 0; jPart < path.nPart; jPart += 1) {
        path.bead(jPart,iBead) -> storeNodeDistance();  // Store old nodal distances      
        N0 += path.getN(jPart,iBead);  // Calculate old nodal action
      }
    }
  }  

  // Perform the bisection (move exactly through kinetic action)
  int skip;
  double tauEff, sigma, sigma2;
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

      VA[iLevel] += path.getV(beadB)*skip;

      rng.normRand(dr, 0, sigma);
      beadB -> r = 0.5 * (beadA -> r + beadC -> r) + dr;
      //path.PutInBox(beadB -> r);

      VB[iLevel] += path.getV(beadB)*skip; 
          
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

  // Check constraint
  if(path.fermi) {
    for (unsigned int iBead = nodeBead0; iBead < nodeBead1; iBead += 1) {
      if(!path.checkConstraint(iBead)) {
        // Restore coordinates
        path.restoreR(affBeads);
        return 0;     
      }    
    }
  }

  double N1(0.0);
  if(path.useNodeDist) {
    // New nodal action
    for (unsigned int iBead = nodeBead0; iBead < nodeBead1; iBead += 1) {
      for (unsigned int jPart = 0; jPart < path.nPart; jPart += 1) {
        path.updateNodeDistance(jPart,iBead);  // Update nodal distances 
        N1 += path.getN(jPart,iBead);  // Calculate new nodal action      
      }
    } 
    
    // Sample nodal action
    if ((N0 - N1) < log(rng.unifRand()))  {
    
      // Restore coordinates
      path.restoreR(affBeads); 
           
      // Restore nodal distances
      for (unsigned int iBead = nodeBead0; iBead < nodeBead1; iBead += 1) {
        for (unsigned int jPart = 0; jPart < path.nPart; jPart += 1) {
          path.bead(jPart,iBead) -> restoreNodeDistance();          
        }
      }
      return 0;    
    }     
  }
   
  // Move Accepted
  return 1; 
}
