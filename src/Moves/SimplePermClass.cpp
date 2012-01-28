#include "SimplePermClass.h"

void SimplePerm::MakeMove()
{
  for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) { 
    nAccept += DoSimplePerm(iBead);
    nAttempt++;
  }
}

// Make a permutation move of 3 particles 
int SimplePerm::DoSimplePerm( const int iBead )
{
  if(path.nPart < 3) return 0; // Only set up for 3 particle permutations
  if((path.fermi)&&(!iBead)) return 0; // Can't permute reference bead

  setPermRadius(iBead); // Set permutation radius

  path.allAction(path.storeRp); // Store original coordinates
  path.allAction(path.storePartRecordp); // Store original particle permutation record

  double K0 = 0.0;
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) 
    K0 += path.getK(iPart, iBead);  // Calculate original kinetic action

  double A0 = K0; // Total affected original action

  // Make a permutation move (try all possible ones in timeslice)
  int permType;
  double K1, A1;
  for (unsigned int i = 0; i < path.nPart - 2; i += 1) {
    for (unsigned int j = i + 1; j < path.nPart - 1; j += 1) {
      for (unsigned int k = j + 1; k < path.nPart; k += 1) {
        if(checkPermRadius(i,j,k)) {
          permType = permute(iBead,i,j,k,-1); // -1 indicates random permutation
           
          K1 = 0.0;
          for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) 
            K1 += path.getK(iPart, iBead);  // Calculate new kinetic action
          
          A1 = K1; // Total affected new action
  
          // Decide whether or not to accept move
          if ((A0 - A1) > log(rng.unifRand()))  {
            path.permCount(permType,1) += 1;
            assignParticleLabels();
            return 1; 
          } else  {
            path.allAction(path.restoreRp); // Restore original coordinates
            path.allAction(path.restorePartRecordp); // Restore original particle permutation record
          }
        }
      }
    }
  }
  return 0;  
}

// Permute Paths (3 atoms at a time)
int SimplePerm::permute( const int slice , const int i , const int j , const int k , int permType )
{
  int perm[path.nPart], iPerm[path.nPart];  

  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    perm[iPart] = iPart;
    iPerm[iPart] = iPart;
  }
  
  if (permType < 0) {
    if (!path.fermi) permType = rng.unifRand(6)-1;  // Bosons     
    else permType = rng.unifRand(3)-1;  // Fermions
  }

  path.permCount(permType,0) += 1;   
  
  setPerm(permType,perm,iPerm,i,j,k);
    
  // Switch bead values according to chosen permutation, i
  for (unsigned int iBead = slice; iBead < path.nBead; iBead += 1)  {
    path.bead(i,iBead) -> r = path.bead(perm[i],iBead) -> rC; 
    path.bead(j,iBead) -> r = path.bead(perm[j],iBead) -> rC; 
    path.bead(k,iBead) -> r = path.bead(perm[k],iBead) -> rC; 
  }
  //std::cout << i << " " << j << " " << k << " " << permType << "\n";
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {  
    path.bead(iPart,path.nBead-1) -> next = path.bead( perm[iPart] , path.nBead-1 ) -> nextC; 
    path.bead(iPart,0) -> prev = path.bead( iPerm[path.bead(iPart,0) -> prevC -> p] , path.nBead-1 ); 
  }
  
  return permType;
}

// Set permutation radius
void SimplePerm::setPermRadius( const int iBead )
{
  double sep2, permRadius2 = path.lambda*iBead*path.tau;

  for (unsigned int iPart = 0; iPart < path.nPart - 1; iPart += 1)  {
    for (unsigned int jPart = iPart + 1; jPart < path.nPart; jPart += 1) {
      dr = path.bead(iPart,iBead) -> r - path.bead(jPart,iBead) -> r;
      path.PutInBox(dr);
      sep2 = dot( dr , dr );
      if(sep2 < permRadius2) {
        path.seps(iPart,jPart) = 1;             
        path.seps(jPart,iPart) = 1;           
      } else {
        path.seps(iPart,jPart) = 0;
        path.seps(jPart,iPart) = 0;
      }
    }
  }
}

// Check permutation radius
bool SimplePerm::checkPermRadius( const int i , const int j , const int k )
{
  if(path.seps(i,j)&&path.seps(i,k)&&path.seps(j,k)) return 1;
  return 0;
}
