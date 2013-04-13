#include "DisplaceBeadClass.h"

void DisplaceBead::MakeMove()
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      nAccept += DoDisplaceBead(iPart,iBead);
      nAttempt++;
    }
  }
}

// Make a single bead displace move
int DisplaceBead::DoDisplaceBead( const int iPart , const int iBead )
{
  Bead *b = path.bead(iPart,iBead);

  // Store coordinates
  b -> storeR();

  path.mode = 0;
  double V0 = path.getV(b);  // Calculate old potential action 
  double K0 = path.getK(b);  // Calculate old kinetic action
  double A0 = V0 + K0;   // Total old action

  // Make a bead move
  rng.unifRand(dr, stepSize);
  b -> move(dr);

  path.mode = 1;
  double V1 = path.getV(b);  // Calculate new potential action
  double K1 = path.getK(b);  // Calculate new kinetic action
  double A1 = V1 + K1;   // Total new action

  // Decide whether or not to accept move
  if ((A0 - A1) < log(rng.unifRand())) {
    b -> restoreR(); // Restore coordinates
    return 0;
  }

  return 1;
}

/*
// Make a single bead displace move
int DisplaceBead::DoDisplace( const int iPart , const int iBead )
{ 
  std::vector<Bead*> beads;
  beads.push_back(path -> bead[iPart][iBead]);

  double V0 = 0.0;  // old potential action 
  double V1 = 0.0;  // new potential action
  double K0 = 0.0;  // old kinetic action
  double K1 = 0.0;  // new kinetic action

  for (std::vector<Bead*>::const_iterator bead = beads.begin(); bead != beads.end(); ++bead) {
    
    // Store coordinates
    (*bead) -> storeR();

    // Calculate old action
    V0 += path -> getV((*bead));
    K0 += path -> getK((*bead));
  
    // Make a bead move
    for (unsigned int iD = 0; iD < path -> nD; iD += 1) dr[iD] = normRand(0,1);
    dr.normalize();
    dr *= stepSize;  
    (*bead) -> move(dr);

    // Calculate new action
    V1 += path -> getV((*bead));
    K1 += path -> getK((*bead));
  }  

  double A0 = V0 + K0;  // Total old action 
  double A1 = V1 + K1;   // Total new action  
  
  // Decide whether or not to accept move  
  if ((A0 - A1) < log(unifRand())) {
    for (std::vector<Bead*>::const_iterator bead = beads.begin(); bead != beads.end(); ++bead) {
      // Restore coordinates
      (*bead) -> restoreR();
    }
    return 0;    
  }
  
  return 1; 
}*/
