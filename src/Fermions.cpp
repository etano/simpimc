#include "PathClass.h"

// Store Rho
void Path::storeRho( const int iBead )
{
  rhoC.slice(iBead) = rho.slice(iBead);
}

// Restore Rho
void Path::restoreRho( const int iBead )
{
  rho.slice(iBead) = rhoC.slice(iBead);
}

// Update Nodal Distances
void Path::updateNodeDistance( const int iPart , const int iBead )
{
  double rmag, r0mag, diff;

  if (!iPart) {
    vec dr = (bead(iPart,iBead) -> r - bead(iPart+1,iBead) -> r);
    bead(iPart,iBead) -> nDist = std::abs(norm(dr,2));
  } else {
    vec dr = (bead(iPart,iBead) -> r - bead(iPart-1,iBead) -> r);
    bead(iPart,iBead) -> nDist = std::abs(norm(dr,2));
  }

  //if(!iBead) bead(iPart,iBead) -> nDist = 0.0; // <<<<<<--------------------Verify
  //else {
  //  gradRho = rho.slice(iBead);
  //  switch (nodeType){
  //    case 1:
  //      for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
  //        dr = (bead(iPart,iBead) -> r) - (bead(jPart,0) -> r);
  //        PutInBox(dr);
  //        diff = norm( dr , 2 );
  //        gradRho(iPart,jPart) = cf2(1,iBead) * diff * rho(iPart,jPart,iBead);
  //      }
  //      break;
  //    default:
  //      for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
  //        rmag = norm( bead(iPart,iBead) -> r , 2 );
  //        r0mag = norm( bead(jPart,0) -> r , 2 );
  //        gradRho(iPart,jPart) = -2.0 * cf2(0,iBead) * (r0mag - rmag*cf3(0,iBead)) * rho(iPart,jPart,iBead);
  //      }
  //      break;
  //  }
  //  bead(iPart,iBead) -> nDist = std::abs( det(rho.slice(iBead)) / det(gradRho) );
  //}
}

// Update Nodal Distances
void Path::updateNodeDistance( std::vector<Bead*>& beads )
{
  for (std::vector<Bead*>::const_iterator b = beads.begin(); b != beads.end(); ++b)
    updateNodeDistance( (*b) -> p , (*b) -> b );
}

// Update Rho
void Path::updateRho ( const int iBead )
{
  double sum, prod, diff;

  switch (nodeType){
    case 1:
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
          dr = (bead(iPart,iBead) -> r) - (bead(jPart,0) -> r);
          PutInBox(dr);
          diff = dot( dr , dr );
          rho(iPart,jPart,iBead) = exp( cf1(1,iBead) * diff );
        }
      }
      break;
    default:
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
          sum = dot( bead(iPart,0)->r , bead(iPart,0)->r ) + dot( bead(jPart,iBead)->r , bead(jPart,iBead)->r );
          prod = dot( bead(iPart,0)->r , bead(jPart,iBead)->r );
          rho(iPart,jPart,iBead) = cf1(0,iBead)*exp(cf2(0,iBead)*(cf3(0,iBead)*sum - 2.*prod));
        }
      }
      break;
  }

}

bool Path::checkConstraint( const int iBead )
{
  int sign = 1.0; // <<<<<<--------------------FIX

  if (!iBead) { // Reference bead moves
    for (unsigned int jBead = 1; jBead < nBead; jBead += 1) { 
      updateRho(jBead); // Update rho
      detRho(jBead) = det(rho.slice(jBead));
      if ( detRho(jBead)*sign < 0.0 ) return 0;
    }
  } else { // Non-reference bead moves
    updateRho(iBead); // Update rho
    detRho(iBead) = det(rho.slice(iBead));
    if ( detRho(iBead)*sign < 0.0 ) return 0;
  }

  return 1;
  //if (bead(0,iBead) -> r[0] >= bead(1,iBead) -> r[0])
  //  return 0;
  //else
  //  return 1;
}

// Check the Constraint 
bool Path::checkConstraint( std::vector<int>& slices )
{
  int sign = 1.0; // <<<<<<--------------------FIX

  for (std::vector<int>::const_iterator iBead = slices.begin(); iBead != slices.end(); ++iBead) {
    updateRho(*iBead); // Update rho
    detRho(*iBead) = det(rho.slice(*iBead));
    if ( detRho(*iBead)*sign < 0.0 ) return 0;
  }

  return 1;
}
