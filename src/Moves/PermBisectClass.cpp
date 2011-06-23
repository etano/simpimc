#include "PermBisectClass.h"

void PermBisect::MakeMove()
{
  nLevel = int(stepSize); // Level of bisection
  nBisectBeads = pow(2,nLevel); // Number of beads in bisection

  nAccept += DoPermBisect();
  nAttempt++;
}

// Make a bisecting permutation move of 3 particles 
int PermBisect::DoPermBisect()
{
  if(path.nPart < 3) return 0; // Can only do 3 particle permuations

  unsigned int bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  unsigned int bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle
  double permTot0 = constructPermTable(bead0,bead1,nBisectBeads,rollOver); // Permutation weight table
    
  int permParts[3], permType;
  permType = selectPerm(permParts,permTot0); // Select Permutation 
  
  // Set up pointers
  Bead *beadA[3], *beadB[3], *beadC[3], *beadI[3], *beadF[3];
  for (unsigned int i = 0; i < 3; i += 1) {
    beadI[i] = path.bead(permParts[i],bead0);
    beadF[i] = beadI[i] -> nextB(nBisectBeads);
  }
  
  // Old nodal action
  double N0(0.0), N1(0.0);
  if(path.useNodeDist) {
    if(rollOver) {
      N0 = 0.0;
      for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
        for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
          path.bead(iPart,iBead) -> storeNodeDistance();  // Store old nodal distances      
          N0 += path.getN(iPart,iBead);  // Calculate old nodal action
        }
      }
    } else {
      for (unsigned int iBead = bead0; iBead < bead1 + 1; iBead += 1) {
        for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
          path.bead(iPart,iBead) -> storeNodeDistance();  // Store old nodal distances      
          N0 += path.getN(iPart,iBead);  // Calculate old nodal action
        }
      }    
    }
  }
  
  // Store permutation record 
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    path.bead(iPart,path.bL(bead1)) -> storePartRecord();    
    path.bead(iPart,path.bL(bead1-1)) -> storePartRecord();      
  }
  
  // Permute particles  
  permuteb(beadF,permType);  
  
  // Set final bisection beads and store
  affBeads.clear();
  for (unsigned int i = 0; i < 3; i += 1) {
    beadF[i] = beadI[i];
    for (unsigned int j = 0; j < nBisectBeads; j += 1) {
      affBeads.push_back(beadF[i]);
      beadF[i] -> storeR();
      beadF[i] = beadF[i] -> next;
    }
  }
  
  // Perform bisection. Move exactly through free action
  int skip;
  double tauEff, sigma;
  double VA[nLevel], VB[nLevel], dA[nLevel+1];
  double dAold = 0.0;  
  dA[nLevel] = 0.0; 
  for (int iLevel = nLevel-1; iLevel >= 0; iLevel -= 1) { 
    skip = pow(2,iLevel);
    tauEff = path.tau*skip;
    sigma = sqrt(lambda*tauEff);
    VA[iLevel] = 0.0;
    VB[iLevel] = 0.0;
    
    for (unsigned int i = 0; i < 3; i += 1) {
      beadA[i] = beadI[i];   
      while(beadA[i] != beadF[i]) {
        beadB[i] = beadA[i] -> nextB(skip); 
        beadC[i] = beadB[i] -> nextB(skip);
              
        VA[iLevel] += path.getV(beadB[i])*skip;

        rng.normRand(dr, 0, sigma);
        beadB[i] -> r = 0.5*(beadA[i] -> r + beadC[i] -> r) + dr;
//        path.PutInBox(beadB[i]);

        VB[iLevel] += path.getV(beadB[i])*skip;
        
        beadA[i] = beadC[i];
      }    
    }

    dA[iLevel] = VA[iLevel] - VB[iLevel];
    dAold = 0.5*(dAold + dA[iLevel+1]);
    // Decide whether or not to accept move (note: only need potential action)
    if ((dA[iLevel] - dAold) < log(rng.unifRand()))  {
    
      // Restore original coordinates
      path.restoreR(affBeads);
      
      // Restore permutation record
      for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
        path.bead(iPart,path.bL(bead1)) -> restorePartRecord();    
        path.bead(iPart,path.bL(bead1-1)) -> restorePartRecord();         
      }   

      return 0;  
    } 
  }
  
  // Construct Permutation Table 
  double permTot1 = constructPermTable(bead0,bead1,nBisectBeads,rollOver);  
 
  // Decide whether or not to accept whole bisection
  if ((permTot0/permTot1) < rng.unifRand())  {
  
    // Restore original coordinates
    path.restoreR(affBeads);
      
    // Restore permutation record
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
      path.bead(iPart,path.bL(bead1)) -> restorePartRecord();    
      path.bead(iPart,path.bL(bead1-1)) -> restorePartRecord();   
    }
    
    return 0;
  }
  
  // Sample constraint and Nodal action
  if(path.fermi) {
    // Check constraint
    if(rollOver) {
      for (unsigned int iBead = 1; iBead < path.nBead; iBead += 1) {
        if(!path.checkConstraint(iBead)) {
    
          // Restore original coordinates
          path.restoreR(affBeads);
          
          // Restore permutation record
          for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
            path.bead(iPart,path.bL(bead1)) -> restorePartRecord();    
            path.bead(iPart,path.bL(bead1-1)) -> restorePartRecord();      
          } 
          
          return 0;     
        }  
      }
    } else { 
      for(beadA[0] = beadI[0]; beadA[0] != beadF[0]; beadA[0] = beadA[0] -> next) {
        if(!path.checkConstraint(beadA[0] -> b)) {
    
          // Restore original coordinates
          path.restoreR(affBeads);
          
          // Restore permutation record
          for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
            path.bead(iPart,path.bL(bead1)) -> restorePartRecord();    
            path.bead(iPart,path.bL(bead1-1)) -> restorePartRecord();     
          }  
          
          return 0;     
        }    
      }
    }
  }
   
  if(path.useNodeDist) { 
    // New nodal action
    N1 = 0.0;
    if(rollOver) {
      for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
        for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
          path.updateNodeDistance(iPart,iBead);  // Update nodal distances 
          N1 += path.getN(iPart,iBead);  // Calculate new nodal action  
        }
      }
    } else {
      for (unsigned int iBead = bead0; iBead < bead1 + 1; iBead += 1) {
        for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
          path.updateNodeDistance(iPart,iBead);  // Update nodal distances    
          N1 += path.getN(iPart,iBead);  // Calculate old nodal action
        }
      }     
    }

    // Sample nodal action
    if ((N0 - N1) < log(rng.unifRand()))  { // Decide whether or not to accept move
    
      // Restore original coordinates
      path.restoreR(affBeads);
      
      // Restore permutation record
      for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
        path.bead(iPart,path.bL(bead1)) -> restorePartRecord();    
        path.bead(iPart,path.bL(bead1-1)) -> restorePartRecord();   
      }
          
      // Restore particle labels
      assignParticleLabels();  

      // Restore nodal distances 
      if(rollOver) {
        for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
          for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
            path.bead(iPart,iBead) -> restoreNodeDistance();          
          }
        }
      } else {
        for (unsigned int iBead = bead0; iBead < bead1 + 1; iBead += 1) {
          for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
            path.bead(iPart,iBead) -> restoreNodeDistance();   
          }
        }      
      }
      
      return 0;    
    }  
  }

  // Assign particle labels
  // HACK!!!
  assignParticleLabels(); // HACK!!!
  // HACK!!!

  path.permCount(permType,1) += 1;  
  
  return 1;    
}

// Construct Permutation Table for 3 Particle Exchanges
double PermBisect::constructPermTable( const int bead0 , const int bead1 , const int nBisectBeads , const bool rollOver ) 
{ 
  int perm[path.nPart], iPerm[path.nPart];
  double cofactor = path.oneOver4LamTau/(1.0*nBisectBeads);
  double permTot = 0.0;
  double diff;
  
  Bead *b0[path.nPart], *b1[path.nPart];
  int b1ip, b1jp, b1kp, n = 0;
  
  if(rollOver) {
    for (unsigned int i = 0; i < path.nPart; i += 1) {
      b0[i] = path.bead(i,bead0);
      b1[i] = b0[i] -> nextB(nBisectBeads);
    } 
  } else {
    for (unsigned int i = 0; i < path.nPart; i += 1) {
      b0[i] = path.bead(i,bead0);
      b1[i] = path.bead(i,bead1);
    } 
  }
  
  for (unsigned int permType = 0; permType < path.nPermType; permType += 1) {
    for (unsigned int i = 0; i < path.nPart - 2; i += 1) {
      for (unsigned int j = i + 1; j < path.nPart - 1; j += 1) {
        for (unsigned int k = j + 1; k < path.nPart; k += 1) { 
        
          b1ip = b1[i] -> p;
          b1jp = b1[j] -> p;
          b1kp = b1[k] -> p;
        
          // Set permutation
          setPerm(permType,perm,iPerm,b1ip,b1jp,b1kp);
          
          // Calculate weight
          diff = 0.0;
          dr = b0[i] -> r - b1[perm[b1ip]] -> r;
          diff += dot( dr , dr );
          dr = b0[j] -> r - b1[perm[b1jp]] -> r;
          diff += dot( dr , dr );
          dr = b0[k] -> r - b1[perm[b1kp]] -> r;
          diff += dot( dr , dr );
          
          permTable(n) = exp(-diff*cofactor);
          permTot += permTable(n);
          n += 1;
        }
      }
    }
  }
  
  return permTot;
}

int PermBisect::selectPerm( int* permParts , double permTot )
{
  double permSubTot = 0.0;
  double x = rng.unifRand(0.0,permTot);  
  
  int n = 0;
  for (unsigned int permType = 0; permType < path.nPermType; permType += 1) {
    for (unsigned int i = 0; i < path.nPart - 2; i += 1) {
      for (unsigned int j = i + 1; j < path.nPart - 1; j += 1) {
        for (unsigned int k = j + 1; k < path.nPart; k += 1) {   
          permSubTot += permTable(n);
          if (x < permSubTot) {
            permParts[0] = i;
            permParts[1] = j;
            permParts[2] = k;
            return permType;
          }          
          n += 1;
        }
      }
    }
  }  
  
  // In case something goes wrong  
  std::cout << "selectPerm Messed Up!" << " " << permTot << " " << x << " ";
  
  return 0;
}

// Permute Paths (3 atoms at a time)
unsigned int PermBisect::permuteb( Bead *b[3] , int permType )
{  
  // permuation holders
  int perm[3], iPerm[3];  
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    perm[iPart] = iPart;
    iPerm[iPart] = iPart;
  }
  
  // Get permutation type
  if (permType < 0) {
    if (path.fermi) permType = rng.unifRand(3)-1;  // Fermions 
    else permType = rng.unifRand(6)-1;  // Bosons    
  }
  path.permCount(permType,0) += 1;   
  
  // Set permutation type
  setPerm(permType,perm,iPerm,0,1,2);
  
  // Assign permutation
  Bead *Pb[3], *iPb[3];  
  for (unsigned int i = 0; i < 3; i += 1) {
    Pb[i] = b[perm[i]];
    iPb[i] = b[iPerm[i]];
  }
  
  // Execute the permutation
  for (unsigned int i = 0; i < 3; i += 1) {
    b[i] -> prev -> next = Pb[i];
    b[i] -> prev = iPb[i] -> prevC;
  }  
    
  return permType;
}
