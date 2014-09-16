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
  vol = pow(L,nD);

  // Constants
  tau = beta/(1.*nBead);

  // Initialize Particles
  nPart = 0;
  vector<Input> speciesInput = in.getChild("Particles").getChildList("Species");
  nSpecies = speciesInput.size();
  for (int iS=0; iS<nSpecies; iS+=1) {
    speciesList.push_back(new Species(speciesInput[iS],iS,nBead,nD));
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
    bool allBeads = in.getChild("Init").getAttribute<bool>("allBeads",false);
    fstream initFile;
    initFile.open(fileName.c_str(), std::ios_base::in);
    for (int iP=0; iP<nPart; ++iP) {
      if (allBeads) {
        for (int iB=0; iB<nBead; ++iB) {
          initFile >> iP >> iB;
          for (int iD=0; iD<nD; ++iD)
            initFile >> bead(iP,iB)->r(iD);
          bead(iP,iB)->storeR();
        }
      } else {
        initFile >> iP;
        Tvector r(nD);
        for (int iD=0; iD<nD; ++iD)
          initFile >> r(iD);
        for (int iB=0; iB<nBead; ++iB) {
          bead(iP,iB)->r = r;
          bead(iP,iB)->storeR();
        }
      }
    }
    initFile.close();
  } else if (initType == "Random") {
    for (int iP=0; iP<nPart; ++iP) {
      for (int iD=0; iD<nD; ++iD) {
        for (int iB=0; iB<nBead; ++iB) {
          RealType tmpRand = rng.unifRand();
          bead(iP,iB)->r(iD) = tmpRand;
        }
      }
      for (int iB=0; iB<nBead; ++iB)
        bead(iP,iB)->storeR();
    }
  }
   // for (int iP=0; iP<nPart; ++iP) {
   //   for (int iB=0; iB<nBead; ++iB) {
   //     cout << iP << " " << iB << " ";
   //     for (int iD=0; iD<nD; ++iD) {
   //        cout << bead(iP,iB)->r(iD) << " ";
   //     }
   //     cout << endl;
   //   }
   // }

   // Reset kCut
   kC = 0.;

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

// Avoid negative and 0 k vectors
inline bool Path::Include(Tvector &k, RealType kCut)
{
  RealType k2 = dot(k,k);
  if (k2 < kCut*kCut && k2 != 0.) {
    if(k(0) > 0.)
      return true;
    else if (nD >= 2 && (k(0) == 0. && k(1) > 0.))
      return true;
    else if (nD >= 3 && (k(0) == 0. && k(1) == 0. && k(2) > 0.))
      return true;
    else
      return false;
  } else
    return false;
}

// Setup universal k vectors according to k cutoff value
// Note that this can be made more general with combinations + permutations
void Path::SetupKs(RealType kCut)
{
  if (kCut <= kC)
    return;
  else {
    ks.clear();
    kIndices.clear();
    magKs.clear();
  }

  // Calculate k box
  kBox.set_size(nD);
  for (int iD=0; iD<nD; iD++)
    kBox(iD) = 2.*M_PI/L;

  // Calculate max k index based on k cutoff
  maxKIndex.set_size(nD);
  for (int iD=0; iD<nD; iD++)
    maxKIndex(iD) = (int) ceil(1.1*kCut/kBox(iD));

  // Set up C vector
  C.resize(nD);
  for (int iD=0; iD<nD; iD++)
    C[iD].set_size(2*maxKIndex(iD)+1);

  // Generate all possible combinations and permutations of k indices
  vector<int> indices;
  for (int iD=0; iD<nD; iD++)
    for (int i=-maxKIndex(iD); i<=maxKIndex(iD); i++)
      indices.push_back(i);
  vector< vector<int> > tmpKIndices;
  genCombPermK(tmpKIndices, indices, nD, false, true);

  // Iterate through indices, form k vectors, and include only those that should be included
  for (int iK=0; iK<tmpKIndices.size(); iK++) {
    Tvector k(nD);
    Ivector ki(nD);
    for (int iD=0; iD<nD; iD++) {
      k(iD) = tmpKIndices[iK][iD]*kBox(iD);
      ki(iD) = maxKIndex(iD) + tmpKIndices[iK][iD];
    }
    if (Include(k, kCut)) {
      ks.push_back(k);
      kIndices.push_back(ki);
      magKs.push_back(sqrt(dot(k,k)));
    }
  }

  cout << "kCut: " << kCut << ", # k vecs: " << ks.size() << endl;

  // Resize rhoK
  rhoK.set_size(nBead,nSpecies,ks.size());
  rhoKC.set_size(nBead,nSpecies,ks.size());
  UpdateRhoK();

}

inline Ccube& Path::GetRhoK()
{
  if (mode)
    return (rhoK);
  else
    return (rhoKC);
}

void Path::CalcRhoKs(int iB, string species)
{
  // Retrieve old or new rho k
  Ccube &tmpRhoK = GetRhoK();

  // Get species info
  int iS, offset;
  GetSpeciesInfo(species,iS,offset);

  // Zero out rho_k array
  for (int iK=0; iK<kIndices.size(); iK++)
    tmpRhoK(beadLoop[iB],iS,iK) = 0.;

  // Calculate rho_k values
  for (int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
    Tvector r = GetR((*this)(iP,iB));
    for (int iD=0; iD<nD; iD++) {
      ComplexType tmpC;
      RealType phi = r(iD)*kBox(iD);
      tmpC = ComplexType(cos(phi), sin(phi));
      C[iD](maxKIndex(iD)) = 1.;
      for (int iK=1; iK<=maxKIndex(iD); iK++) {
        C[iD](maxKIndex(iD)+iK) = tmpC * C[iD](maxKIndex[iD]+iK-1);
        C[iD](maxKIndex(iD)-iK) = conj(C[iD](maxKIndex[iD]+iK));
      }
    }
    for (int iK=0; iK<kIndices.size(); iK++) {
      ComplexType factor = 1;
      for (int iD=0; iD<nD; iD++)
        factor *= C[iD](kIndices[iK](iD));
      tmpRhoK(beadLoop[iB],iS,iK) += factor;
    }
  }
}

// Update all rho k values
void Path::UpdateRhoK()
{
  bool tmpMode = GetMode();
  SetMode(0);
  for (int iB=0; iB<nBead; iB++) {
    for (int iS=0; iS<nSpecies; iS++)
      CalcRhoKs(iB, speciesList[iS]->name);
  }
  rhoK = rhoKC;

  SetMode(tmpMode);
}

// Update rho k values for specific particles and slices
void Path::UpdateRhoK(int b0, int b1, vector<int> &particles, int level)
{
  bool tmpMode = GetMode();
  int skip = 1<<level;

  // Update C values of changed particles
  for (int iP=0; iP<particles.size(); iP++) {
    int p = particles[iP];
    int iS = bead(p,0)->species.iS;
    for (int iB=b0; iB<=b1; iB+=skip) {
      // Reset new to old
      for (int iK=0; iK<kIndices.size(); iK++)
        rhoK(beadLoop[iB],iS,iK) = rhoKC(beadLoop[iB],iS,iK);

      // Add in new values
      SetMode(1);
      Tvector rNew = GetR((*this)(p,iB));
      for (int iD=0; iD<nD; iD++) {
        ComplexType tmpC;
        RealType phi = rNew(iD)*kBox(iD);
        tmpC = ComplexType(cos(phi), sin(phi));
        C[iD](maxKIndex(iD)) = 1.;
        for (int iK=1; iK<=maxKIndex(iD); iK++) {
          C[iD](maxKIndex(iD)+iK) = tmpC * C[iD](maxKIndex[iD]+iK-1);
          C[iD](maxKIndex(iD)-iK) = conj(C[iD](maxKIndex[iD]+iK));
        }
      }
      for (int iK=0; iK<kIndices.size(); iK++) {
        Ivector &ki = kIndices[iK];
        ComplexType factor = 1.;
        for (int iD=0; iD<nD; iD++)
          factor *= C[iD](ki(iD));
        rhoK(beadLoop[iB],iS,iK) += factor;
      }

      // Subtract out old values
      SetMode(0);
      Tvector rOld = GetR((*this)(p,iB));
      for (int iD=0; iD<nD; iD++) {
        ComplexType tmpC;
        RealType phi = rOld(iD)*kBox(iD);
        tmpC = ComplexType(cos(phi), sin(phi));
        C[iD](maxKIndex(iD)) = 1.;
        for (int iK=1; iK<=maxKIndex(iD); iK++) {
          C[iD](maxKIndex(iD)+iK) = tmpC * C[iD](maxKIndex[iD]+iK-1);
          C[iD](maxKIndex(iD)-iK) = conj(C[iD](maxKIndex[iD]+iK));
        }
      }
      for (int iK=0; iK<kIndices.size(); iK++) {
        Ivector &ki = kIndices[iK];
        ComplexType factor = 1.;
        for (int iD=0; iD<nD; iD++)
          factor *= C[iD](ki(iD));
        rhoK(beadLoop[iB],iS,iK) -= factor;
      }
    }
  }

  SetMode(tmpMode);

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
