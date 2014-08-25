#include "PathClass.h"

void Path::Init(Input &in, IOClass &out, RNG &rng)
{
  nD = in.getChild("System").getAttribute<int>("nD");
  nBead = in.getChild("System").getAttribute<int>("nBead");
  beta = in.getChild("System").getAttribute<RealType>("beta");
  PBC = in.getChild("System").getAttribute<int>("PBC", 1);
  if (PBC)
    L = in.getChild("System").getAttribute<RealType>("L");
  else
    L = 0.;

  // Constants
  tau = beta/(1.*nBead);

  // Initialize Particles
  nPart = 0;
  vector<Input> speciesInput = in.getChild("Particles").getChildList("Species");
  nSpecies = speciesInput.size();
  for (int iS=0; iS<nSpecies; iS+=1) {
    speciesList.push_back(new Species(speciesInput[iS],nBead,nD));
    nPart += speciesList[iS]->nPart;
  }

  // Maximum Bisection Level
  maxLevel = int(log2(nBead))-1;

  // Initiate bead looping
  beadLoop.set_size(2*nBead);
  for (unsigned int iB = 0; iB < nBead; iB += 1) {
    beadLoop(iB) = iB;
    beadLoop(iB + nBead) = beadLoop(iB);
  }

  // Initiate mode
  SetMode(1);

  // Initiate beads
  bead.set_size(nPart,nBead);
  for (unsigned int iS=0; iS<nSpecies; iS+=1)
    for (unsigned int iP=0; iP<speciesList[iS]->nPart; iP+=1)
      for (unsigned int iB=0; iB<nBead; iB+=1)
        bead(iS*speciesList[iS]->nPart + iP,iB) = new Bead(nD,*speciesList[iS],iS*speciesList[iS]->nPart + iP,iB);

  // Initiate bead connections
  for (unsigned int iP = 0; iP < nPart; iP += 1) {
    bead(iP,0) -> next = bead(iP,1);
    bead(iP,0) -> prev = bead(iP,nBead-1);
    for (unsigned int iB = 1; iB < nBead; iB += 1) {
      bead(iP,iB) -> next = bead(iP,beadLoop(iB+1));
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

  // Initialize paths
  string initType = in.getChild("Init").getAttribute<string>("type");
  if (initType == "File") {
    string fileName = in.getChild("Init").getAttribute<string>("name");
    fstream initFile;
    initFile.open(fileName.c_str(), std::ios_base::in);
    for (int iB=0; iB<nBead; ++iB) {
      for (int iP=0; iP<nPart; ++iP) {
        initFile >> iB >> iP;
        for (int iD=0; iD<nD; ++iD)
          initFile >> bead(iP,iB)->r(iD);
        bead(iP,iB) -> storeR();
      }
    }
    initFile.close();
  } else if (initType == "Random") {
    for (int iP=0; iP<nPart; ++iP) {
      for (int iD=0; iD<nD; ++iD) {
        RealType tmpRand = rng.unifRand();
        for (int iB=0; iB<nBead; ++iB)
          bead(iP,iB)->r(iD) = tmpRand;
      }
      for (int iB=0; iB<nBead; ++iB)
        bead(iP,iB) -> storeR();
    }
  }

}

// Get species info
void Path::GetSpeciesInfo(string species, int &iSpecies, int &offset)
{
  int tmpOffset = 0;
  for (unsigned int iS=0; iS<nSpecies; iS+=1) {
    if (speciesList[iS]->name == species) {
      iSpecies = iS;
      offset = tmpOffset;
    }
    tmpOffset += speciesList[iS]->nPart;
  }
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
