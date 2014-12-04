#include "PathClass.h"

void Path::Init(Input &in, IOClass &out, RNG &rng)
{
  out.CreateGroup("System");
  nD = in.getChild("System").getAttribute<int>("nD");
  out.Write("System/nD",nD);
  nBead = in.getChild("System").getAttribute<int>("nBead");
  out.Write("System/nBead",nBead);
  beta = in.getChild("System").getAttribute<double>("beta");
  out.Write("System/beta",beta);
  PBC = in.getChild("System").getAttribute<int>("PBC", 1);
  out.Write("System/PBC",PBC);
  if (PBC) {
    L = in.getChild("System").getAttribute<double>("L");
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

  // Initialize species
  out.CreateGroup("System/Particles");
  vector<Input> speciesInput = in.getChild("Particles").getChildList("Species");
  nSpecies = speciesInput.size();
  nPart = 0;
  for (int iS=0; iS<nSpecies; iS++) {
    speciesList.push_back(new Species(speciesInput[iS],out,iS,nD,nBead));
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

  // Initialize paths within each species
  for (int iS=0; iS<nSpecies; iS++)
    speciesList[iS]->InitPaths(speciesInput[iS],out,rng,InterComm,L);

  // Reset kCut
  kC = 0.;
  double kCut = in.getChild("System").getAttribute<double>("kCut",2.*M_PI/(pow(vol,1./nD)));
  SetupKs(kCut);

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
  for (int iS=0; iS<nSpecies; ++iS) {
    for (int iP=0; iP<speciesList[iS]->nPart; ++iP) {
      for (int iB=0; iB<nBead; ++iB) {
        cout << iP << " " << iB << " ";
        vec<double> r = GetR((*this)(iS,iP,iB));
        for (int iD=0; iD<nD; ++iD) {
           cout << r(iD) << " ";
        }
        cout << endl;
      }
    }
  }
  cout << endl;
}

// Get species info
void Path::GetSpeciesInfo(string species, int &iSpecies)
{
  for (unsigned int iS=0; iS<nSpecies; iS++) {
    if (speciesList[iS]->name == species) {
      iSpecies = iS;
      return;
    }
  }
  cout << "ERROR: No species of name " << species << " !" << endl;
}

// Put R in the Box
void Path::PutInBox(vec<double>& r)
{
  for (int iD=0; iD<nD; ++iD)
    r(iD) -= nearbyint(r(iD)*iL)*L;
}

// Avoid negative and 0 k vectors
bool Path::Include(vec<double> &k, double kCut)
{
  double k2 = dot(k,k);
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
void Path::SetupKs(double kCut)
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
    vec<double> k(nD);
    vec<int> ki(nD);
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

  //cout << "kCut: " << kCut << ", # k vecs: " << ks.size() << endl;

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
      for (int iP=0; iP<speciesList[iS]->nPart; ++iP) {
        (*this)(iS,iP,iB)->rhoK.zeros(kIndices.size());
        (*this)(iS,iP,iB)->rhoKC.zeros(kIndices.size());
      }
    }
  }

  // Calculate rho_k's
  for (int iB=0; iB<nBead; iB++) {
    for (int iS=0; iS<nSpecies; iS++) {
      // Calculate rho_kp values
      for (int iP=0; iP<speciesList[iS]->nPart; iP++) {
        CalcRhoKP((*this)(iS,iP,iB));
        (*this)(iS,iP,iB) -> restoreRhoK();
        rhoKC(iB,iS) += (*this)(iS,iP,iB)->rhoKC;
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

      // Calculate rho_kp values
      for (int iP=0; iP<speciesList[iS]->nPart; iP++) {
        rhoK(iB,iS) += (*this)(iS,iP,iB)->rhoK;
      }
    }
  }

}

// Update rho k values for specific particles and slices
void Path::UpdateRhoKP(int b0, int b1, vector< pair<int,int> > &particles, int level)
{
  bool tmpMode = GetMode();
  SetMode(1);
  int skip = 1<<level;

  // Get species numbers
  vector<int> species, ps;
  for (int p=0; p<particles.size(); p++) {
    int iS = particles[p].first;
    int iP = particles[p].second;
    species.push_back(iS);
    ps.push_back(iP);
  }

  for (int s=0; s<species.size(); ++s)
    UpdateRhoKP(b0, b1, species[s], ps, level);

  SetMode(tmpMode);

}

// Update rho k values for specific particles and slices
void Path::UpdateRhoKP(int b0, int b1, int iS, vector<int> &particles, int level)
{
  bool tmpMode = GetMode();
  SetMode(1);
  int skip = 1<<level;

  // Reset to old copy
  for (int iB=b0; iB<b1; iB+=skip)
    rhoK(beadLoop(iB),iS) = rhoKC(beadLoop(iB),iS);

  // Update C values of changed particles
  for (int p=0; p<particles.size(); p++) {
    int iP = particles[p];
    for (int iB=b0; iB<b1; iB+=skip) {
      // Calculate new values
      CalcRhoKP((*this)(iS,iP,iB));

      // Add in new values and subtract out old values
      rhoK(beadLoop(iB),iS) += (*this)(iS,iP,iB)->rhoK - (*this)(iS,iP,iB)->rhoKC;

    }
  }

  SetMode(tmpMode);

}

// Update rho k values for specific particles and slices
void Path::UpdateRhoK(int b0, int b1, vector< pair<int,int > > &particles, int level)
{
  bool tmpMode = GetMode();
  int skip = 1<<level;

  // Get species numbers
  vector<int> species;
  for (int p=0; p<particles.size(); p++) {
    int iS = particles[p].first;
    int iP = particles[p].second;
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
    int iS = particles[p].first;
    int iP = particles[p].second;
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
void Path::CalcC(vec<double> &r)
{
  for (int iD=0; iD<nD; iD++) {
    complex<double> tmpC;
    double phi = r(iD)*kBox(iD);
    tmpC = complex<double>(cos(phi), sin(phi));
    C(iD)(maxKIndex(iD)) = 1.;
    for (int iK=1; iK<=maxKIndex(iD); iK++) {
      C(iD)(maxKIndex(iD)+iK) = tmpC * C(iD)(maxKIndex(iD)+iK-1);
      C(iD)(maxKIndex(iD)-iK) = conj(C(iD)(maxKIndex(iD)+iK));
    }
  }
}

// Add rho_k for a single particle
void Path::AddRhoKP(field< vec< complex<double> > >& tmpRhoK, int iP, int iB, int iS, int pm)
{
  vec<double> r = GetR((*this)(iS,iP,iB));
  CalcC(r);
  for (int iK=0; iK<kIndices.size(); iK++) {
    vec<int> &ki = kIndices[iK];
    complex<double> factor = pm;
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
  vec<double> r = GetR(b);
  vec< complex<double> >& tmpRhoK = GetRhoK(b);
  CalcC(r);
  for (int iK=0; iK<kIndices.size(); iK++) {
    complex<double> factor = 1.;
    for (int iD=0; iD<nD; iD++)
      factor *= C(iD)(kIndices[iK](iD));
    tmpRhoK(iK) = factor;
  }
}

// Compute path sign
int Path::CalcSign()
{
  sign = 1;
  for (int iS=0; iS<nSpecies; iS++) {
    if (speciesList[iS]->fermi) {
      vector<int> cycles;
      SetCycleCount(iS, cycles);
      for (int i=0; i<cycles.size(); i++)
        sign *= pow(-1,cycles[i]-1);
    }
  }
}

// Retrieve list of cycles
void Path::SetCycleCount(int iS, vector<int>& cycles)
{
  GetSpeciesInfo(speciesList[iS]->name,iS);
  vec<int> alreadyCounted;
  alreadyCounted.zeros(speciesList[iS]->nPart);
  for (unsigned int iP=0; iP<speciesList[iS]->nPart; iP++) {
    if (!alreadyCounted(iP)) {
      int cycleLength = 1;
      Bead* b = (*this)(iS,iP,nBead-1);
      alreadyCounted(iP) = 1;
      while (GetNextBead(b,1) != (*this)(iS,iP,0)) {
        cycleLength++;
        b = GetNextBead(b,nBead);
        alreadyCounted(b->p) = 1;
      }
      cycles.push_back(cycleLength);
    }
  }
}

// Retrieve permutation sector
int Path::GetPermSector(int iS)
{
  vector<int> cycles;
  SetCycleCount(iS, cycles);
  return GetPermSector(iS, cycles);
}

// Retrieve permutation sector
int Path::GetPermSector(int iS, vector<int>& cycles)
{
  sort(cycles.begin(),cycles.end());
  possPermsIterator = possPerms.find(cycles);
  if (possPermsIterator == possPerms.end()) {
    cout << "Broken Permutation: " << endl;
    for (vector<int>::size_type i=0; i != cycles.size(); i++)
      cout << cycles[i] << " ";
    cout << endl;
    exit(1);
  }
  return possPermsIterator->second;
}

// Setup permutation sectors
void Path::SetupPermSectors(int n, int sectorsMax)
{
  cout << "Setting up permutation sectors..." << endl;
  vector<int> a;
  a.resize(n);
  for (int i=0; i<n; i++) {
    a[i] = 0;
  }
  int k = 1;
  int y = n-1;
  vector< vector<int> > tmpPossPerms;
  while (k != 0 && (sectorsMax > possPerms.size() || !sectorsMax)) {
    int x = a[k-1] + 1;
    k -= 1;
    while (2*x <= y) {
      a[k] = x;
      y -= x;
      k += 1;
    }
    int l = k+1;
    while (x <= y && (sectorsMax > possPerms.size() || !sectorsMax)) {
      a[k] = x;
      a[l] = y;
      vector<int> b;
      for (vector<int>::size_type j=0; j!=k+2; j++)
        b.push_back(a[j]);
      tmpPossPerms.push_back(b);
      x += 1;
      y -= 1;
    }
    a[k] = x+y;
    y = x+y-1;
    vector<int> c;
    for (vector<int>::size_type j=0; j!=k+1; j++)
      c.push_back(a[j]);
    tmpPossPerms.push_back(c);
  }

  int nSectors = tmpPossPerms.size();
  for (vector<int>::size_type j=0; j != nSectors; j++)
    possPerms[tmpPossPerms[j]] = j;

}
