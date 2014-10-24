#include "PathClass.h"

void Path::Init(Input &in, IOClass &out, RNG &rng)
{
  out.CreateGroup("System");
  nD = in.getChild("System").getAttribute<int>("nD");
  out.Write("System/nD",nD);
  nBead = in.getChild("System").getAttribute<int>("nBead");
  out.Write("System/nBead",nBead);
  beta = in.getChild("System").getAttribute<RealType>("beta");
  out.Write("System/beta",beta);
  PBC = in.getChild("System").getAttribute<int>("PBC", 1);
  out.Write("System/PBC",PBC);
  if (PBC) {
    L = in.getChild("System").getAttribute<RealType>("L");
    iL = 1./L;
    vol = pow(L,nD);
  } else {
    L = 0.;
    iL = 0.;
    vol = 1.;
  }
  out.Write("System/L",L);

  // Approximate with fast math
  approximate = false;

  // Constants
  tau = beta/(1.*nBead);
  out.Write("System/tau",tau);

  // Initialize Particles
  out.CreateGroup("System/Particles");
  vector<Input> speciesInput = in.getChild("Particles").getChildList("Species");
  nSpecies = speciesInput.size();
  nPart = 0;
  for (int iS=0; iS<nSpecies; iS++) {
    speciesList.push_back(new Species(speciesInput[iS],out,iS,nBead,nD));
    nPart += speciesList[iS]->nPart;
  }
  out.Write("System/nPart",nPart);

  // Maximum Bisection Level
  maxLevel = int(log2(nBead))-1;

  // Initiate bead looping
  beadLoop.set_size(2*nBead);
  for (unsigned int iB = 0; iB < nBead; iB++) {
    beadLoop(iB) = iB;
    beadLoop(iB + nBead) = beadLoop(iB);
  }

  // Initiate mode
  SetMode(1);

  // Initiate beads
  bead.set_size(nPart,nBead);
  unsigned int unique_id = 0;
  for (int iS=0; iS<nSpecies; iS++) {
    int offset;
    GetSpeciesInfo(speciesList[iS]->name,iS,offset);
    for (unsigned int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
      for (unsigned int iB=0; iB<nBead; iB++) {
        bead(iP,iB) = new Bead(nD,*speciesList[iS],iP,iB,unique_id);
        unique_id++;
      }
    }
  }

  // Initiate bead connections
  for (unsigned int iP = 0; iP < nPart; iP++) {
    bead(iP,0) -> next = bead(iP,1);
    bead(iP,0) -> nextC = bead(iP,1);
    bead(iP,0) -> prev = bead(iP,nBead-1);
    bead(iP,0) -> prevC = bead(iP,nBead-1);
    for (unsigned int iB = 1; iB < nBead; iB++) {
      bead(iP,iB) -> next = bead(iP,beadLoop(iB+1));
      bead(iP,iB) -> nextC = bead(iP,beadLoop(iB+1));
      bead(iP,iB) -> prev = bead(iP,iB-1);
      bead(iP,iB) -> prevC = bead(iP,iB-1);
    }
  }


  // Set permutation types
  iCount.zeros(nPart);
  pCount.zeros(nPart);
  setPType();

  // Initialize paths
  InitPaths(in,out,rng);

  // Reset kCut
  kC = 0.;
  if (PBC) {
    RealType kCut = in.getChild("System").getAttribute<RealType>("kCut");
    SetupKs(kCut);
  }

  // Initiate nodal things
  sign = 1;
  refBead = 0;
  CalcSign();
}

// Store R
void Path::storeR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
    (*beadIter) -> storeR();
}

// Restore R
void Path::restoreR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
    (*beadIter) -> restoreR();
}

void Path::PrintPath()
{
  SetMode(0);
  for (int iP=0; iP<nPart; ++iP) {
    for (int iB=0; iB<nBead; ++iB) {
      cout << iP << " " << iB << " ";
      Tvector r = GetR(bead(iP,iB));
      for (int iD=0; iD<nD; ++iD) {
         cout << r(iD) << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
}

// Get species info
void Path::GetSpeciesInfo(string species, int &iSpecies, int &offset)
{
  int tmpOffset = 0;
  for (unsigned int iS=0; iS<nSpecies; iS++) {
    if (speciesList[iS]->name == species) {
      iSpecies = iS;
      offset = tmpOffset;
      break;
    }
    tmpOffset += speciesList[iS]->nPart;
  }
}

// Put R in the Box
void Path::PutInBox(Tvector& r)
{
  Tvector n(NDIM);
  n(0) = -nearbyint(r(0)*iL);
  r(0) += n(0)*L;
#if NDIM>1
  n(1) = -nearbyint(r(1)*iL);
  r(1) += n(1)*L;
#endif
#if NDIM>2
  n(2) = -nearbyint(r(2)*iL);
  r(2) += n(2)*L;
#endif
}

// Avoid negative and 0 k vectors
bool Path::Include(Tvector &k, RealType kCut)
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
void Path::SetupKs(RealType kCut)
{
  if (kCut <= kC)
    return;
  else {
    ks.clear();
    kIndices.clear();
    magKs.clear();
    kC = kCut;
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
  C.set_size(nD);
  for (int iD=0; iD<nD; iD++)
    C(iD).set_size(2*maxKIndex(iD)+1);

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

  // Initiate rho k
  InitRhoK();

}

// Update all rho k values
void Path::InitRhoK()
{

  // Resize rhoK
  rhoK.set_size(nBead,nSpecies);
  rhoKC.set_size(nBead,nSpecies);

  // Set to old mode
  bool tmpMode = GetMode();
  SetMode(0);

  // Zero out rho_k array
  for (int iB=0; iB<nBead; iB++) {
    for (int iS=0; iS<nSpecies; iS++) {
      rhoK(iB,iS).zeros(kIndices.size());
      rhoKC(iB,iS).zeros(kIndices.size());
      for (int iP=0; iP<nPart; ++iP) {
        bead(iP,iB)->rhoK.zeros(kIndices.size());
        bead(iP,iB)->rhoKC.zeros(kIndices.size());
      }
    }
  }

  // Calculate rho_k's
  for (int iB=0; iB<nBead; iB++) {
    for (int iS=0; iS<nSpecies; iS++) {
      // Get species info
      int offset;
      GetSpeciesInfo(speciesList[iS]->name,iS,offset);

      // Calculate rho_kp values
      for (int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
        CalcRhoKP(bead(iP,iB));
        bead(iP,iB)->restoreRhoK();
        rhoKC(iB,iS) += bead(iP,iB)->rhoKC;
      }
      rhoK(iB,iS) = rhoKC(iB,iS);
    }
  }

  SetMode(tmpMode);
}

// Update all rho k values
void Path::UpdateRhoK()
{
  // Calculate rho_k's
  for (int iB=0; iB<nBead; iB++) {
    for (int iS=0; iS<nSpecies; iS++) {
      // Zero out rho_k array
      rhoK(iB,iS).zeros();

      // Get species info
      int offset;
      GetSpeciesInfo(speciesList[iS]->name,iS,offset);

      // Calculate rho_kp values
      for (int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
        rhoK(iB,iS) += bead(iP,iB)->rhoK;
      }
    }
  }

}

// Update rho k values for specific particles and slices
void Path::UpdateRhoKP(int b0, int b1, vector<int> &particles, int level)
{
  bool tmpMode = GetMode();
  SetMode(1);
  int skip = 1<<level;

  // Get species numbers
  vector<int> species;
  for (int p=0; p<particles.size(); p++) {
    int iP = particles[p];
    int iS = bead(iP,beadLoop(b0))->species.iS;
    species.push_back(iS);
  }

  // Reset to old copy
  for (int s=0; s<species.size(); s++) {
    int iS = species[s];
    for (int iB=b0; iB<b1; iB+=skip)
      rhoK(beadLoop(iB),iS) = rhoKC(beadLoop(iB),iS);
  }

  // Update C values of changed particles
  for (int p=0; p<particles.size(); p++) {
    int iP = particles[p];
    int iS = bead(iP,0)->species.iS;
    for (int iB=b0; iB<b1; iB+=skip) {
      // Calculate new values
      SetMode(1);
      CalcRhoKP(bead(iP,beadLoop(iB)));

      // Add in new values and subtract out old values
      rhoK(beadLoop(iB),iS) += bead(iP,beadLoop(iB))->rhoK - bead(iP,beadLoop(iB))->rhoKC;

    }
  }

  SetMode(tmpMode);

}

// Update rho k values for specific particles and slices
void Path::UpdateRhoK(int b0, int b1, vector<int> &particles, int level)
{
  bool tmpMode = GetMode();
  int skip = 1<<level;

  // Get species numbers
  vector<int> species;
  for (int p=0; p<particles.size(); p++) {
    int iP = particles[p];
    int iS = bead(iP,beadLoop(b0))->species.iS;
    species.push_back(iS);
  }

  // Reset to old copy
  for (int s=0; s<species.size(); s++) {
    int iS = species[s];
    for (int iB=b0; iB<b1; iB+=skip)
      rhoK(beadLoop(iB),iS) = rhoKC(beadLoop(iB),iS);
  }

  // Update C values of changed particles
  for (int p=0; p<particles.size(); p++) {
    int iP = particles[p];
    int iS = bead(iP,0)->species.iS;
    for (int iB=b0; iB<b1; iB+=skip) {
      // Add in new values
      SetMode(1);
      AddRhoKP(rhoK,iP,iB,iS,1);

      // Subtract out old values
      SetMode(0);
      AddRhoKP(rhoK,iP,iB,iS,-1);
    }
  }

  SetMode(tmpMode);

}

// Calculate C values for rho_k
void Path::CalcC(Tvector &r)
{
  for (int iD=0; iD<nD; iD++) {
    ComplexType tmpC;
    RealType phi = r(iD)*kBox(iD);
    tmpC = ComplexType(cos(phi), sin(phi));
    C(iD)(maxKIndex(iD)) = 1.;
    for (int iK=1; iK<=maxKIndex(iD); iK++) {
      C(iD)(maxKIndex(iD)+iK) = tmpC * C(iD)(maxKIndex(iD)+iK-1);
      C(iD)(maxKIndex(iD)-iK) = conj(C(iD)(maxKIndex(iD)+iK));
    }
  }
}

// Add rho_k for a single particle
void Path::AddRhoKP(field<Cvector>& tmpRhoK, int iP, int iB, int iS, int pm)
{
  Tvector r = GetR((*this)(iP,beadLoop(iB)));
  CalcC(r);
  for (int iK=0; iK<kIndices.size(); iK++) {
    Ivector &ki = kIndices[iK];
    ComplexType factor = pm;
    for (int iD=0; iD<nD; iD++)
      factor *= C(iD)(ki(iD));
    tmpRhoK(beadLoop(iB),iS)(iK) += factor;
  }
}

// Store rhoK
void Path::storeRhoKP(vector<Bead*>& affBeads)
{
  if (kC > 0)
    for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
      (*beadIter) -> storeRhoK();
}

// Restore rhoK
void Path::restoreRhoKP(vector<Bead*>& affBeads)
{
  if (kC > 0)
    for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
      (*beadIter) -> restoreRhoK();
}

// Calc rho_k for a single particle
inline void Path::CalcRhoKP(Bead* b)
{
  Tvector r = GetR(b);
  Cvector& tmpRhoK = GetRhoK(b);
  CalcC(r);
  for (int iK=0; iK<kIndices.size(); iK++) {
    Ivector &ki = kIndices[iK];
    ComplexType factor = 1.;
    for (int iD=0; iD<nD; iD++)
      factor *= C(iD)(ki(iD));
    tmpRhoK(iK) = factor;
  }
}

// Compute path sign
int Path::CalcSign()
{
  // Count up each permutation type
  sign = 1;
  for (int iS=0; iS<nSpecies; iS++) {
    if (speciesList[iS]->fermi) {
      int offset;
      GetSpeciesInfo(speciesList[iS]->name,iS,offset);
      Ivector pCount;
      pCount.zeros(speciesList[iS]->nPart);
      for (unsigned int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
        iCount(iP) = 0;
        Bead* b = bead(iP,nBead-1);
        while (b->next != bead(iP,0)) {
          iCount(iP)++;
          b = bead(b->next->p,nBead-1);
        }
        pCount(iCount(iP))++; // Bin values by number of exchanges
      }
      for (unsigned int iP=0; iP<speciesList[iS]->nPart; iP++) {
        if (iP%2)
          sign *= pow(-1,pCount(iP));
      }
    }
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
