#include "PathClass.h"

Path::Path( const int nPartIn, const int nDIn , const int nBeadIn, const double betaIn , const double lambdaIn , const int fermiIn , const int halfspaceIn , const int nodeTypeIn , const int useNodeDistIn , const double LIn )
  : nPart(nPartIn) , nD(nDIn) , nBead(nBeadIn) , beta(betaIn) , lambda(lambdaIn) , fermi(fermiIn) , halfspace(halfspaceIn) , nodeType(nodeTypeIn) , useNodeDist(useNodeDistIn) , L(LIn) , oneOverL(1.0/LIn)
{
  // Constants
  tau = beta/(1.0*nBead);
  oneOverLamTau = 1.0/(1.0*lambda*tau);
  oneOver4LamTau = 1.0/(4.0*lambda*tau);
  oneOver4LamTau2 = 1.0/(4.0*lambda*tau*tau);
  nPartnBeadnDOver2Tau = (1.0*nPart*nBead*nD)/(2.0*tau);
  halfTauOmega2 = 0.5*tau*omega*omega;
  halfOmega2 = 0.5*omega*omega;
    
  // Initiate bead looping
  bL.set_size(2*nBead);  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    bL(iBead) = iBead;
    bL(iBead + nBead) = bL(iBead);
  }
  
  // Initiate beads
  bead.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      bead(iPart,iBead) = new Bead(nD,iPart,iBead);
    }
  }
       
  // Initiate bead connections
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    bead(iPart,0) -> next = bead(iPart,1);
    bead(iPart,0) -> prev = bead(iPart,nBead-1);
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      bead(iPart,iBead) -> next = bead(iPart,bL(iBead+1));
      bead(iPart,iBead) -> prev = bead(iPart,iBead-1);    
    }
  }

  // Set permutation types
  if (fermi) nPermType = 3;
  else nPermType = 6; 
  iCount.zeros(nPart);
  pCount.zeros(nPart); 
  setPType();

  // Initiate permutation counter
  permCount.zeros(nPermType,2);

  // Initiate separation matrix
  seps.zeros(nPart,nPart);

  // Set up cofactors for nodal calculations
  cf1.set_size(3,nBead);
  cf2.set_size(3,nBead);

  switch (nodeType){
    case 1:
      for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
        cf1(1,iBead) = -1.0/(4.0*lambda*iBead*tau);
        cf2(1,iBead) = -2.0/(4.0*lambda*iBead*tau);
      }
      break;
    default:
      cf1(0,0) = 0;
      cf2(0,0) = 0;
      for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
        cf1(0,iBead) = omega/sinh(iBead*tau*omega);
        cf2(0,iBead) = cosh(iBead*tau*omega);
      }
      break;
  }

  // Initiate rho, rhoC
  rho.zeros(nPart,nPart,nBead);
  rhoC.zeros(nPart,nPart,nBead);
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    updateRho(iBead);
    storeRho(iBead);
  }

  // Initiate detRho
  detRho.zeros(nBead);

  // Initiate gradRho
  gradRho.zeros(nPart,nPart);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      bead(iPart,iBead) -> nDist = 1.50;
    }
  }  
  
  // Bead Member Functions
  storeRp = &Bead::storeR;
  restoreRp = &Bead::restoreR;
  storePartRecordp = &Bead::storePartRecord;
  restorePartRecordp = &Bead::restorePartRecord;
  storeNodeDistancep = &Bead::storeNodeDistance;
  restoreNodeDistancep = &Bead::restoreNodeDistance;

  // Error Counts
  infCount = 0;
  nanCount = 0;
  errCount = 0;
}

// Identify permuation state
int Path::getPType()
{
  if (nType < 2) return 0;

  unsigned int iType, xPart;

  // Count up each permutation type
  pCount.zeros();
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {   
    iCount(iPart) = 0;  
    xPart = iPart;
    while ((bead(xPart,nBead-1) -> next -> p) != iPart) {
      iCount(iPart) += 1;
      xPart = (bead(xPart,nBead-1) -> next -> p);
    }
    pCount(iCount(iPart)) += 1; // Bin values by number of exchanges
  }    
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) 
    pCount(iPart) = pCount(iPart)/(iPart+1); // Normalize values

  iType = 0;
  while (iType < nType) {
    xPart = 0;
    while (xPart < nPart) {
      if (pCount(xPart) != pType(iType,xPart)) break;
      xPart += 1;
    }
    if (xPart == nPart) break;
    iType += 1;      
  }

  return iType;
}

// Set up permuation states
void Path::setPType()
{
  switch (nPart) {
    case 1:
      nType = 1;
      pType.zeros(nType,nPart);
      pType(0,0) = 1; 
      break;
    case 2:
      nType = 2;
      pType.zeros(nType,nPart);
      pType(0,0) = 2;
      pType(1,1) = 1;
      break;
    case 3:
      nType = 3;
      pType.zeros(nType,nPart);
      pType(0,0) = 3;
      pType(1,2) = 1; 
      pType(2,0) = 1;
      pType(2,1) = 1;
      break;
    case 4:
      nType = 5;
      pType.zeros(nType,nPart);
      pType(0,0) = 4;
      pType(1,0) = 2;
      pType(1,1) = 1;
      pType(2,0) = 1;
      pType(2,2) = 1;
      pType(3,1) = 2;
      pType(4,3) = 1;
      break;
    case 5:
      nType = 7;
      pType.zeros(nType,nPart);
      pType(0,0) = 5;
      pType(1,0) = 3;
      pType(1,1) = 1;
      pType(2,0) = 2;
      pType(2,2) = 1;
      pType(3,0) = 1;
      pType(3,3) = 1;
      pType(4,0) = 1;
      pType(4,1) = 2;
      pType(5,1) = 1;
      pType(5,2) = 1;
      pType(6,4) = 1;
      break;
    default:
      nType = 1;
      pType.zeros(nType,nPart);
      pType(0,0) = 1;
      break;
  }
}

// Print out permutation configuration
void Path::printPerm()
{  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1)
    std::cout << bead(iPart,nBead-1) -> self -> p << " " << bead(iPart,nBead-1) -> next -> p << "\n";
  std::cout << "\n";    
}

// Print out bead positions
void Path::printBeads ()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      std::cout << iPart << " " << iBead << " " << bead(iPart,iBead) -> r << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

