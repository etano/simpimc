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

  if (nPart == 1) return;

  if (nD == 1) {
    double dr(0.0), dr1(0.0), dr2(0.0);
    if (mode) {
      if (iPart==0) {
        dr = (bead(1,iBead) -> r[0] - bead(0,iBead) -> r[0]);
      } else if(iPart==nPart-1) {
        dr = (bead(nPart-1,iBead) -> r[0] - bead(nPart-2,iBead) -> r[0]);
      } else {
        dr1 = (bead(iPart+1,iBead) -> r[0] - bead(iPart,iBead) -> r[0]);
        dr2 = (bead(iPart,iBead) -> r[0] - bead(iPart-1,iBead) -> r[0]);
        dr = std::min(dr1,dr2);
      }
      bead(iPart,iBead) -> nDist = dr;
    } else {
      if (iPart==0) {
        dr = (bead(1,iBead) -> rC[0] - bead(0,iBead) -> rC[0]);
      } else if(iPart==nPart-1) {
        dr = (bead(nPart-1,iBead) -> rC[0] - bead(nPart-2,iBead) -> rC[0]);
      } else {
        dr1 = (bead(iPart+1,iBead) -> rC[0] - bead(iPart,iBead) -> rC[0]);
        dr2 = (bead(iPart,iBead) -> rC[0] - bead(iPart-1,iBead) -> rC[0]);
        dr = std::min(dr1,dr2);
      }
      bead(iPart,iBead) -> nDistC = dr;
    }
  } else {

    if(!iBead) bead(iPart,iBead) -> nDist = 0.0; // <<<<<<--------------------Verify
    else {
      gradRho = rho.slice(iBead);
      if (mode) {
        switch (nodeType){
          case 1:
            for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
              dr = (bead(iPart,iBead) -> r) - (bead(jPart,0) -> r);
              PutInBox(dr);
              diff = norm( dr , 2 );
              gradRho(iPart,jPart) = cf2(1,iBead) * diff * rho(iPart,jPart,iBead);
            }
            break;
          default:
            for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
              rmag = norm( bead(iPart,iBead) -> r , 2 );
              r0mag = norm( bead(jPart,0) -> r , 2 );
              gradRho(iPart,jPart) = -2.0 * cf2(0,iBead) * (r0mag - rmag*cf3(0,iBead)) * rho(iPart,jPart,iBead);
            }
            break;
        }
        bead(iPart,iBead) -> nDist = std::abs( det(rho.slice(iBead)) / det(gradRho) );
      } else {
        switch (nodeType){
          case 1:
            for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
              dr = (bead(iPart,iBead) -> rC) - (bead(jPart,0) -> rC);
              PutInBox(dr);
              diff = norm( dr , 2 );
              gradRho(iPart,jPart) = cf2(1,iBead) * diff * rho(iPart,jPart,iBead);
            }
            break;
          default:
            for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
              rmag = norm( bead(iPart,iBead) -> rC , 2 );
              r0mag = norm( bead(jPart,0) -> rC , 2 );
              gradRho(iPart,jPart) = -2.0 * cf2(0,iBead) * (r0mag - rmag*cf3(0,iBead)) * rho(iPart,jPart,iBead);
            }
            break;
        }
        bead(iPart,iBead) -> nDistC = std::abs( det(rho.slice(iBead)) / det(gradRho) );
      }
    }
  }
}

// Update Nodal Distances
void Path::updateNodeDistance( Bead *b )
{
  updateNodeDistance( b-> p , b -> b );
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
  int sliceDiff = int(std::min(double(iBead),double(nBead-iBead)));

  switch (nodeType){
    case 1:
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
          dr = (bead(iPart,iBead) -> rC) - (bead(jPart,0) -> rC);
          PutInBox(dr);
          diff = dot( dr , dr );
          rho(iPart,jPart,iBead) = exp( cf1(1,sliceDiff) * diff );
        }
      }
      break;
    default:
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int jPart = 0; jPart < nPart; jPart += 1) {
          sum = dot( bead(iPart,0)->rC , bead(iPart,0)->rC ) + dot( bead(jPart,iBead)->rC , bead(jPart,iBead)->rC );
          prod = dot( bead(iPart,0)->rC , bead(jPart,iBead)->rC );
          rho(iPart,jPart,iBead) = cf1(0,sliceDiff)*exp(cf2(0,sliceDiff)*(cf3(0,sliceDiff)*sum - 2.*prod));
        }
      }
      break;
  }

}

bool Path::checkConstraint( const int iBead )
{
  if (nPart == 1) return 1;
  if (nD == 1) {
    if (mode) {
      for (unsigned int iPart = 1; iPart < nPart; iPart += 1) {
        if (bead(iPart,iBead) -> r[0] < bead(iPart-1,iBead) -> r[0])
          return 0;
      }
      return 1;
    } else {
       for (unsigned int iPart = 1; iPart < nPart; iPart += 1) {
        if (bead(iPart,iBead) -> rC[0] < bead(iPart-1,iBead) -> rC[0])
          return 0;
      }
      return 1;
    }
  } else {
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
  }
  return 1;

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
