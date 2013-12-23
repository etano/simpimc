#ifndef PathClass_H
#define PathClass_H

#include "config.h"       // Standard libraries
#include "SpeciesClass.h"
#include "BeadClass.h"
#include "IO/InputClass.h"
#include "IO/IOClass.h"

typedef void (Bead::*BeadMemFn)();

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

class Path
{
private:

protected:

public:
  // Constructor
  Path() {};
  void Init(Input &in, IOClass &out);

  // Given Global Constants
  unsigned int nPart, nD, nBead;
  RealType beta, L, iL;

  // Calculated Global Constants
  RealType tau;
  unsigned int nPermType;
  unsigned int maxLevel;

  // Species
  unsigned int nSpecies;
  vector<Species*> speciesList;

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
  void Dr(Bead* b0, Bead* b1, Tvector &dr);

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(Tvector& r);

  // Mode (use copy or true)
  bool mode;
};

inline Bead* Path::operator() (int iP, int iB)
{
  return bead(iP,beadLoop[iB]);
}

// Put R in the Box
inline void Path::PutInBox(Tvector& r)
{
  if(!PBC) {
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
inline void Path::Dr(Bead* b0, Bead* b1, Tvector &dr)
{
  dr = b0->r - b1->r;
  PutInBox(dr);
}

#endif
