#ifndef PathClass_H
#define PathClass_H

#include "config.h"       // Standard libraries
#include "SpeciesClass.h"
#include "BeadClass.h"
#include "IO/InputFile.h"
#include "IO/IO.h"

typedef void (Bead::*BeadMemFn)();

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

class Path
{
private:
  // Rho Variables
  arma::cube rho, rhoC;
  Tmatrix gradRho;
  Tvector detRho;
  Tmatrix cf1, cf2, cf3;
protected:
  // protected things
public:
  // Constructor
  Path() {};
  void Init(Input &in, IOClass &out);

  // Given Global Constants
  unsigned int nPart, nD, nBead;
  RealType beta;
  RealType L, iL;

  // Species
  unsigned int nSpecies;
  vector<Species*> speciesList;

  // Trapping potential
  bool trap;
  RealType omega;

  // Fermions
  int halfspace;
  int nodeType;
  bool useNodeDist;

  // Calculated Global Constants
  RealType tau;
  unsigned int nPermType;
  RealType nPartnBeadnDOver2Tau, halfTauOmega2, halfOmega2, onePlusTau2Omega2Over12, onePlus3Tau2Omega2Over12;
  unsigned int maxLevel;

  // Permutation Counter
  int getPType();
  void setPType();
  Imatrix permCount;
  Imatrix pType;
  Ivector iCount, pCount;
  unsigned int nType;

  // Print things
  void printPerm();
  void printBeads();

  // Beads
  arma::field<Bead*> bead;

  // Bead Iterator
  std::vector<Bead*>::const_iterator beadIter;

  // Bead Class member functions (redundant)
  inline void storeR( std::vector<Bead*>& affBeads );
  inline void restoreR( std::vector<Bead*>& affBeads );
  BeadMemFn storeRp, restoreRp, storePartRecordp, restorePartRecordp, storeNodeDistancep, restoreNodeDistancep;
  void beadAction( Bead *b , BeadMemFn p );
  void beadsAction( std::vector<Bead*>& beads , BeadMemFn p );
  void partAction( int iP , BeadMemFn p );
  void allAction( BeadMemFn p );

  // Rho Functions
  void storeRho( const int iB );
  void restoreRho( const int iB );
  void updateRho( const int iB );
  void updateNodeDistance( const int iP , const int iB );
  void updateNodeDistance( std::vector<Bead*>& beads );
  void updateNodeDistance( Bead *b );
  bool checkConstraint( const int iB );
  bool checkConstraint( std::vector<int>& slices );

  // Difference Vector
  void Dr(Bead* b0, Bead* b1, Tvector &dr);

  // Periodic Boundary Conditions
  void PutInBox( Tvector& r );

  // Action Functions
  RealType getK();
  RealType getK(const int iP);
  RealType getK(int b0, int b1, vector<int> &particles, int level);
  RealType getV();
  RealType getV( const unsigned int iP );
  RealType getV( const unsigned int iP , const int iB );
  RealType getV( Bead *bi );
  RealType getVint( Bead *b1 , Bead *b2 );
  RealType getVext( Bead *b );
  RealType getN( const int iP , const int iB );
  RealType getN( Bead *b );
  RealType getN( std::vector<Bead*>& beads );
  RealType getN( Bead *b , int skip );
  RealType getN( const int iP , const int iB , int skip );
  RealType getNSlice( const int iB , const int skip );

  // Mode (use copy or true)
  bool mode;

  // Tables
  Ivector bL;
  Tmatrix seps;

  // Bad Number Counts
  unsigned int infCount, nanCount, errCount;
};

// Single bead action
inline void Path::beadAction( Bead *b , BeadMemFn p )
{
  CALL_MEMBER_FN(*b,p)();
}

// Single particle action
inline void Path::partAction( int iP , BeadMemFn p )
{
  for (unsigned int iB = 0; iB < nBead; iB ++)
    CALL_MEMBER_FN(*bead(iP,iB),p)();
}

// Bead set action
inline void Path::beadsAction( std::vector<Bead*>& beads , BeadMemFn p )
{
  for (beadIter = beads.begin(); beadIter != beads.end(); ++beadIter)
    CALL_MEMBER_FN(*(*beadIter),p)();
}

// Whole system action
inline void Path::allAction( BeadMemFn p )
{
  for (unsigned int iP = 0; iP < nPart; iP ++) {
    for (unsigned int iB = 0; iB < nBead; iB ++)
      CALL_MEMBER_FN(*bead(iP,iB),p)();
  }
}

// Put R in the Box
inline void Path::PutInBox( Tvector& r )
{
  if(!trap) {
    for (unsigned int iD = 0; iD < nD; iD++) {
      RealType n = -floor(r(iD) * iL + 0.5);
      r(iD) += n * L;
    }
  }
}

// Store R
inline void Path::storeR( std::vector<Bead*>& affBeads )
{
  for (beadIter = affBeads.begin(); beadIter != affBeads.end(); ++beadIter) {
    (*beadIter) -> storeR();
  }
}

// Restore R
inline void Path::restoreR( std::vector<Bead*>& affBeads )
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
