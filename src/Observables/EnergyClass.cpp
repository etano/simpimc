#include "EnergyClass.h"

void Energy::Accumulate( const int pType )
{
  KE = getKE(); // Kinetic Energy
  VE = getVE(); // Potential Energy  
  NE = getNE(); // Nodal Energy      
  E = KE + VE + NE; // Total Energy
  
  Etot(pType) += E; // Add up total energy
  KEtot(pType) += KE; // Add up total Kinetic energy
  VEtot(pType) += VE; // Add up total Potential energy
  NEtot(pType) += NE; // Add up total Nodal energy     
}

void Energy::Output()
{
  E = 0.0;
  KE = 0.0;
  VE = 0.0;
  NE = 0.0;  
  
  Etot *= oneOverNbeadBlock;
  KEtot *= oneOverNbeadBlock;
  VEtot *= oneOverNbeadBlock;
  NEtot *= oneOverNbeadBlock;
  
  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    trace << Etot(iType) << " " << KEtot(iType) << " " << VEtot(iType) << " " << NEtot(iType) << " ";
    E += Etot(iType);
    KE += KEtot(iType);
    VE += VEtot(iType);
    NE += NEtot(iType);
  }
  trace << E << " " << KE << " " << VE << " " << NE << "\n";

  Etot.zeros();
  KEtot.zeros();
  VEtot.zeros();
  NEtot.zeros();
}

void Energy::Print()
{
  std::cout << "\nE: " << E << "\nKE: " << KE << "\nVE: " << VE << "\nNE: " << NE << "\n"; 
}

void Energy::Stats()
{
  statsEnergy( outputFile.c_str() , path.nType , nBlock );
}

// Get Kinetic Energy
double Energy::getKE()
{
  double tot = 0.0;  
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
      path.dr = path.bead(iPart,iBead) -> r - (path.bead(iPart,iBead) -> next -> r);
      tot += dot( path.dr , path.dr );
    }
  }
  
  return path.nPartnBeadnDOver2Tau - path.oneOver4LamTau2 * tot;
}

// Get Potential Energy
double Energy::getVE()
{
  double VEext = 0.0;
  double VEint = 0.0;  
  for (unsigned int iPart = 0; iPart < path.nPart-1; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
      for (unsigned int jPart = iPart+1; jPart < path.nPart; jPart += 1) {
        VEint += 2.0 * path.getVint( path.bead(iPart,iBead) , path.bead(jPart,iBead) );
      }
      VEext += dot( path.bead(iPart,iBead) -> r , path.bead(iPart,iBead) -> r );
    }
  }
  for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
    VEext += dot( path.bead(path.nPart-1,iBead) -> r , path.bead(path.nPart-1,iBead) -> r );
  }
  VEext *= path.halfOmega2;
  
  return VEext + VEint;
}

// Get Two Bead Interaction Energy
double Energy::getVEint( Bead *b1 , Bead *b2 )
{
  return 0;
}

// Get Nodal Energy
double Energy::getNE()
{
  if(!path.useNodeDist) return 0;

  double xi, tot;
  
  tot = 0.0;
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
      xi = (path.bead(iPart,iBead) -> nDist)*(path.bead(iPart,iBead) -> next -> nDist)*path.oneOverLamTau;
      tot += xi/(path.tau*(exp(xi) - 1.0));
    } 
  }

  return tot;
  
}
