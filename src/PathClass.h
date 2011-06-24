#ifndef PathClass_H
#define PathClass_H

#include "StandardLibraries.h"       // Standard libraries
#include "BeadClass.h"
#include "GlobalConstants.h"

typedef void (Bead::*BeadMemFn)();

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

class Path
{
private:
  // Rho Variables
  cube rho, rhoC;
  mat gradRho;
  vec detRho;
  mat cf1, cf2;
protected:
  // protected things
public:
  // Constructor
  Path( const int nPartIn , const int nDIn , const int nBeadIn, const double betaIn , const int fermiIn , const int halfspaceIn , const int nodeTypeIn , const int useNodeDistIn , const double LIn );  

  // Given Global Constants
  const unsigned int nPart, nD, nBead;
  const double beta;
  const bool fermi;
  const int halfspace;
  const int nodeType;
  const bool useNodeDist;
  const double L;

  // Calculated Global Constants
  double tau;
  unsigned int nPermType;
  double oneOverLamTau, oneOver4LamTau, oneOver4LamTau2, nPartnBeadnDOver2Tau, halfTauOmega2, halfOmega2;
  
  // Permutation Counter
  int getPType();
  void setPType();
  imat permCount;
  imat pType;
  ivec iCount, pCount;
  unsigned int nType;
  
  // Print things
  void printPerm();
  void printBeads();
  
  // Beads
  field<Bead*> bead;

  // Bead Iterator
  std::vector<Bead*>::const_iterator beadIter;

  // Bead Class member functions (redundant)
  inline void storeR( std::vector<Bead*>& affBeads );
  inline void restoreR( std::vector<Bead*>& affBeads );
  BeadMemFn storeRp, restoreRp, storePartRecordp, restorePartRecordp, storeNodeDistancep, restoreNodeDistancep;
  void beadAction( Bead *b , BeadMemFn p );
  void beadsAction( std::vector<Bead*>& beads , BeadMemFn p );
  void partAction( int iPart , BeadMemFn p );
  void allAction( BeadMemFn p );
  
  // Rho Functions
  void storeRho( const int iBead );
  void restoreRho( const int iBead );
  void updateRho( const int iBead );
  void updateNodeDistance( const int iPart , const int iBead );
  void updateNodeDistance( Bead *b );
  void updateNodeDistance( std::vector<Bead*>& beads );
  bool checkConstraint( const int iBead );  
  bool checkConstraint( std::vector<int>& slices );  
  
  // Difference Vector
  vec dr;

  // Periodic Boundary Conditions
  void PutInBox( vec& r );
  
  // Action Functions
  double getK();
  double getK( const int iPart );
  double getK( const int iPart , const int iBead );
  double getK( Bead *bi );
  double getV();
  double getV( const unsigned int iPart );
  double getV( const unsigned int iPart , const int iBead );
  double getV( Bead *bi );
  double getVint( Bead *b1 , Bead *b2 );
  double getN( const int iPart , const int iBead );
  double getN( Bead *b );
  double getN( std::vector<Bead*>& beads );
  
  // Tables
  ivec bL;
  mat seps;
  
  // Bad Number Counts
  unsigned int infCount, nanCount, errCount;
};

// Single bead action
inline void Path::beadAction( Bead *b , BeadMemFn p )
{
  CALL_MEMBER_FN(*b,p)();
}

// Single particle action
inline void Path::partAction( int iPart , BeadMemFn p )
{
  for (unsigned int iBead = 0; iBead < nBead; iBead ++)
    CALL_MEMBER_FN(*bead(iPart,iBead),p)();
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
  for (unsigned int iPart = 0; iPart < nPart; iPart ++) {
    for (unsigned int iBead = 0; iBead < nBead; iBead ++)
      CALL_MEMBER_FN(*bead(iPart,iBead),p)();
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

// Put R in the Box
inline void Path::PutInBox( vec& r )
{
  for (unsigned int iD = 0; iD < nD; iD++) {
    while (r(iD) > L/2) {
      r(iD) -= L;
    }
    while (r(iD) < -L/2) {
      r(iD) += L;
    }
  }
}

#endif
