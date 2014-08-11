#ifndef PathClass_H
#define PathClass_H

#include "SpeciesClass.h"
#include "BeadClass.h"
#include "Utils/config.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"

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

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(Tvector& r);

  // Mode (use copy or true)
  bool mode;
  inline void SetMode(int m) { mode = m; };
};

inline Bead* Path::operator() (int iP, int iB)
{
  return bead(iP,beadLoop[iB]);
}

// Put R in the Box
inline void Path::PutInBox(Tvector& r)
{
  if(PBC) {
    for (unsigned int iD = 0; iD < nD; iD++) {
      RealType n = -floor(r(iD)*iL + 0.5);
      r(iD) += n*L;
    }
  }
}

// Store R
inline void Path::storeR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter) {
    (*beadIter) -> storeR();
  }
}

// Restore R
inline void Path::restoreR(vector<Bead*>& affBeads)
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter) {
    (*beadIter) -> restoreR();
  }
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

#endif
