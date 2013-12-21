#include "PathClass.h"

void Path::Init(Input &in, IOClass &out)
{
  nD = in.get<int>("Input.System.nD");
  nBead = in.get<int>("Input.System.nBead");
  beta = in.get<RealType>("Input.System.beta");
  L = in.get<RealType>("Input.System.L");
  trap = in.get<int>("Input.System.trap", 0);

  //!!!!!!!!!!
  useNodeDist = 0;

  // Initialize Particles
  nPart = 0;
  vector<Input> speciesInput = in.getObjectList("Input.Particles");
  nSpecies = speciesInput.size();
  for (int iS=0; iS<nSpecies; iS+=1) {
    speciesList.push_back(new Species(speciesInput[iS],nBead,nD));
    nPart += speciesList[iS]->nPart;
  }

  // Constants
  tau = beta/(1.*nBead);
  //iLamTau = 1./(1.*speciesList[0]->lambda*tau);
  //i4LamTau = 1./(4.*speciesList[0]->lambda*tau);
  //i4LamTau2 = 1./(4.*speciesList[0]->lambda*tau*tau);
  nPartnBeadnDOver2Tau = (1.*nPart*nBead*nD)/(2.*tau);
  omega = 1.;
  halfTauOmega2 = 0.5*tau*omega*omega;
  halfOmega2 = 0.5*omega*omega;
  onePlusTau2Omega2Over12 = 1. + (tau*tau*omega*omega/12.);
  onePlus3Tau2Omega2Over12 = 1. + (3.*tau*tau*omega*omega/12.);

  // Maximum Bisection Level
  maxLevel = int(log2(nBead))-1;

  // Initiate bead looping
  bL.set_size(2*nBead);
  for (unsigned int iB = 0; iB < nBead; iB += 1) {
    bL(iB) = iB;
    bL(iB + nBead) = bL(iB);
  }

  // Initiate beads
  bead.set_size(nPart,nBead);
  for (unsigned int iS=0; iS<nSpecies; iS+=1) {
    for (unsigned int iP=0; iP<speciesList[iS]->nPart; iP+=1) {
      for (unsigned int iB=0; iB<nBead; iB+=1) {
        bead(iS*speciesList[iS]->nPart + iP,iB) = new Bead(nD,*speciesList[iS],iP,iB);
      }
    }
  }

  // Initiate bead connections
  for (unsigned int iP = 0; iP < nPart; iP += 1) {
    bead(iP,0) -> next = bead(iP,1);
    bead(iP,0) -> prev = bead(iP,nBead-1);
    for (unsigned int iB = 1; iB < nBead; iB += 1) {
      bead(iP,iB) -> next = bead(iP,bL(iB+1));
      bead(iP,iB) -> prev = bead(iP,iB-1);
    }
  }

  // Set permutation types
  if (speciesList[0]->fermi) nPermType = 3;
  else nPermType = 6;
  iCount.zeros(nPart);
  pCount.zeros(nPart);
  setPType();

  // Initiate permutation counter
  permCount.zeros(nPermType,2);

  // Initiate separation matrix
  seps.zeros(nPart,nPart);

  // Set up cofactors for nodal calculations
  cf1.zeros(3,nBead);
  cf2.zeros(3,nBead);
  cf3.zeros(3,nBead);

  switch (0){
    case 1:
      for (unsigned int iB = 1; iB < nBead; iB += 1) {
        cf1(1,iB) = -1.0/(4.0*speciesList[0]->lambda*iB*tau);
        cf2(1,iB) = -2.0/(4.0*speciesList[0]->lambda*iB*tau);
      }
      break;
    default:
      for (unsigned int iB = 1; iB < nBead; iB += 1) {
        cf1(0,iB) = pow((omega/(2.0*pi*sinh(iB*tau*omega))),nD/2.0);
        cf2(0,iB) = -omega/(2.0*sinh(iB*tau*omega));
        cf3(0,iB) = cosh(iB*tau*omega);
      }
      break;
  }

  // Initiate rho, rhoC
  rho.zeros(nPart,nPart,nBead);
  rhoC.zeros(nPart,nPart,nBead);
  for (unsigned int iB = 0; iB < nBead; iB += 1) {
    updateRho(iB);
    storeRho(iB);
  }

  // Initiate detRho
  detRho.zeros(nBead);

  // Initiate gradRho
  gradRho.zeros(nPart,nPart);
  for (unsigned int iP = 0; iP < nPart; iP += 1) {
    for (unsigned int iB = 0; iB < nBead; iB += 1) {
      updateNodeDistance(iP,iB);
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
  for (unsigned int iP = 0; iP < nPart; iP += 1) {
    iCount(iP) = 0;
    xPart = iP;
    while ((bead(xPart,nBead-1) -> next -> p) != iP) {
      iCount(iP) += 1;
      xPart = (bead(xPart,nBead-1) -> next -> p);
    }
    pCount(iCount(iP)) += 1; // Bin values by number of exchanges
  }

  for (unsigned int iP = 0; iP < nPart; iP += 1) 
    pCount(iP) = pCount(iP)/(iP+1); // Normalize values

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
  for (unsigned int iP = 0; iP < nPart; iP += 1)
    std::cout << bead(iP,nBead-1) -> self -> p << " " << bead(iP,nBead-1) -> next -> p << "\n";
  std::cout << "\n";
}

// Print out bead positions
void Path::printBeads ()
{
  for (unsigned int iP = 0; iP < nPart; iP += 1) {
    for (unsigned int iB = 0; iB < nBead; iB += 1) {
      std::cout << iP << " " << iB << " " << bead(iP,iB) -> r << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}
