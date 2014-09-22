#ifndef PathClass_H
#define PathClass_H

#include "SpeciesClass.h"
#include "BeadClass.h"
#include "Utils/config.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"
#include "Utils/Algorithm/Algorithm.h"

class Path
{
private:

protected:

public:
  // Constructor
  Path() {};
  void Init(Input &in, IOClass &out, RNG &rng);

  // Given Global Constants
  unsigned int nPart, nD, nBead;
  RealType beta, L, iL, vol;

  // Calculated Global Constants
  RealType tau;
  unsigned int nPermType;
  unsigned int maxLevel;

  // Species
  unsigned int nSpecies;
  vector<Species*> speciesList;
  void GetSpeciesInfo(string species, int &iSpecies, int &offset);

  // Permutation Counter
  int getPType();
  void setPType();
  Imatrix permCount;
  Imatrix pType;
  Ivector iCount, pCount;
  unsigned int nType;

  // Beads
  arma::field<Bead*> bead;
  vector<Bead*>::const_iterator beadIter;
  Ivector beadLoop;
  inline Bead* operator() (int iP, int iB);
  inline void storeR(vector<Bead*> &affBeads);
  inline void restoreR(vector<Bead*> &affBeads);
  Tvector& GetR(Bead* b);
  void Dr(Tvector &r0, Tvector &r1, Tvector &dr);
  void Dr(Bead* b0, Tvector &r1, Tvector &dr);
  void Dr(Bead* b0, Bead* b1, Tvector &dr);
  void RBar(Bead* b0, Bead* b1, Tvector &rBar);
  void DrDrPDrrP(int b0, int b1, int p0, int p1, RealType &rMag, RealType &rPMag, RealType &rrPMag, Tvector &r0, Tvector &r1, Tvector &dr);
  void PrintPath();

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(Tvector& r);

  // Mode (use copy or true)
  bool mode;
  inline void SetMode(int m) { mode = m; };
  inline bool GetMode() { return mode; };

  // k vectors and rho_k
  vector<Tvector> ks;
  vector<Ivector> kIndices;
  vector<Cvector> C;
  vector<RealType> magKs;
  arma::field<Cvector> rhoK, rhoKC, rhoKP, rhoKPC;
  Tvector kBox;
  RealType kC;
  Ivector maxKIndex;
  void InitRhoK();
  void UpdateRhoK();
  void UpdateRhoK(int b0, int b1, vector<int> &particles, int level);
  void UpdateRhoKP(int b0, int b1, vector<int> &particles, int level);
  void CalcC(Tvector &r);
  void AddRhoKP(arma::field<Cvector>& tmpRhoK, int iP, int iB, int iS, int pm);
  void CalcRhoKP(arma::field<Cvector>& tmpRhoK, int iP, int iB, int iS);
  void SetupKs(RealType kCut);
  inline bool Include(Tvector &k, RealType kCut);
  inline void storeRhoK(Bead* b);
  inline void storeRhoK(vector<Bead*>& affBeads);
  inline void restoreRhoK(vector<Bead*>& affBeads);
  inline void restoreRhoK(Bead* b);
};

inline Bead* Path::operator() (int iP, int iB)
{
  return bead(iP,beadLoop(iB));
}

// Put R in the Box
inline void Path::PutInBox(Tvector& r)
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

// Store R
inline void Path::storeR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
    (*beadIter) -> storeR();
}

// Restore R
inline void Path::restoreR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
    (*beadIter) -> restoreR();
}

// Get dr
inline void Path::RBar(Bead* b0, Bead* b1, Tvector &rBar)
{
  Dr(b0, b1, rBar);
  rBar = GetR(b1) + 0.5*rBar;
}

// Get r
inline Tvector& Path::GetR(Bead* b)
{
  if (mode)
    return (b->r);
  else
    return (b->rC);
}

// Get dr
inline void Path::Dr(Tvector &r0, Tvector &r1, Tvector &dr)
{
  dr = r0 - r1;
  PutInBox(dr);
}

// Get dr, drP, and drrP
inline void Path::DrDrPDrrP(int b0, int b1, int p0, int p1, RealType &rMag, RealType &rPMag, RealType &rrPMag, Tvector &r, Tvector &rP, Tvector &rrP)
{
  b0 = beadLoop(b0);
  b1 = beadLoop(b1);
  r = GetR(bead(p1,b0)) - GetR(bead(p0,b0));
  rP = GetR(bead(p1,b1)) - GetR(bead(p0,b1));

  Tvector n(NDIM), m(NDIM);
  n(0) = -nearbyint(r(0)*iL);
  r(0) += n(0)*L;
  m(0) = nearbyint((r(0)-rP(0))*iL);
  rP(0) += m(0)*L;
#if NDIM>1
  n(1) = -nearbyint(r(1)*iL);
  r(1) += n(1)*L;
  m(1) = nearbyint((r(1)-rP(1))*iL);
  rP(1) += m(1)*L;
#endif
#if NDIM>2
  n(2) = -nearbyint(r(2)*iL);
  r(2) += n(2)*L;
  m(2) = nearbyint((r(2)-rP(2))*iL);
  rP(2) += m(2)*L;
#endif
  rMag = mag(r);
  rPMag = mag(rP);

  rrP = r - rP;
  PutInBox(rrP);
  rrPMag = mag(rrP);
}

// Get dr
inline void Path::Dr(Bead* b0, Tvector &r1, Tvector &dr)
{
  Dr(GetR(b0), r1, dr);
}

// Get dr
inline void Path::Dr(Bead* b0, Bead* b1, Tvector &dr)
{
  Dr(GetR(b0), GetR(b1), dr);
}

// Store rhoK
inline void Path::storeRhoK(Bead* b)
{
  int iB = b->b;
  int iS = b->species.iS;
  int iP = b->p;
  rhoKC(iB,iS) = rhoK(iB,iS);
  rhoKPC(iB,iS,iP) = rhoKP(iB,iS,iP);
}

// Restore rhoK
inline void Path::restoreRhoK(Bead* b)
{
  int iB = b->b;
  int iS = b->species.iS;
  int iP = b->p;
  rhoK(iB,iS) = rhoKC(iB,iS);
  rhoKP(iB,iS,iP) = rhoKPC(iB,iS,iP);
}

// Store rhoK
inline void Path::storeRhoK(vector<Bead*>& affBeads)
{
  if (kC > 0)
    for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
      storeRhoK((*beadIter));
}

// Restore rhoK
inline void Path::restoreRhoK(vector<Bead*>& affBeads)
{
  if (kC > 0)
    for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter)
      restoreRhoK((*beadIter));
}

#endif
