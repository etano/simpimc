#include "PathClass.h"

#include "armadillo"
using namespace arma;

// Get Kinetic Energy
double Path::getKE()
{
  double tot = 0.0;  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)  {
      dr = bead(iPart,iBead) -> r - (bead(iPart,iBead) -> next -> r);
      tot += dot( dr , dr );
    }
  }
  
  return nPartnBeadnDOver2Tau - oneOver4LamTau2*tot;
}

// Get Potential Energy
double Path::getVE()
{
  double VEext = 0.0;
  double VEint = 0.0;  
  for (unsigned int iPart = 0; iPart < nPart-1; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)  {
      for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) {
        VEint += 2.0 * getVint( bead(iPart,iBead) , bead(jPart,iBead) );
      }
      VEext += dot( bead(iPart,iBead) -> r , bead(iPart,iBead) -> r );
    }
  }
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    VEext += dot( bead(nPart-1,iBead) -> r , bead(nPart-1,iBead) -> r );
  }
  VEext *= halfOmega2;
  
  return VEext + VEint;
}

// Get Two Bead Interaction Energy
double Path::getVEint( Bead *b1 , Bead *b2 )
{
  return 0;
}

// Get Nodal Energy
double Path::getNE()
{
  if(!useNodeDist) return 0;

  double xi, tot;
  
  tot = 0.0;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)  {
      xi = (bead(iPart,iBead) -> nDist)*(bead(iPart,iBead) -> next -> nDist)*oneOverLamTau;
      tot += xi/(tau*(exp(xi) - 1.0));
    } 
  }

  return tot;
  
}
