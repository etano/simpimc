#include "EnergyClass.h"

void Energy::Accumulate( const int iType )
{
  KE = getKE(); // Kinetic Energy
  VE = getVE(); // Potential Energy
  NE = getNE(); // Nodal Energy
  E = KE + VE + NE; // Total Energy

  Etot(iType) += E; // Add up total energy
  KEtot(iType) += KE; // Add up total Kinetic energy
  VEtot(iType) += VE; // Add up total Potential energy
  NEtot(iType) += NE; // Add up total Nodal energy
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
  trace << E << " " << KE << " " << VE << " " << NE << endl;

  Etot.zeros();
  KEtot.zeros();
  VEtot.zeros();
  NEtot.zeros();
}

void Energy::Print()
{
  std::cout << "\nE: " << E << "\nKE: " << KE << "\nVE: " << VE << "\nNE: " << NE << endl; 
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
      dr = path.bead(iPart,iBead) -> r - (path.bead(iPart,iBead) -> next -> r);
      path.PutInBox( dr );
      tot += dot( dr , dr );
    }
  }

  return path.nPartnBeadnDOver2Tau - path.oneOver4LamTau2 * tot;
}

// Get Potential Energy
inline double Energy::getVE()
{
  return getVEext() + getVEint();
}

// Get Interaction Potential Energy
double Energy::getVEint()
{
  double tot = 0.0;

  for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
    for (unsigned int iPart = 0; iPart < path.nPart-1; iPart += 1)  {
      for (unsigned int jPart = iPart+1; jPart < path.nPart; jPart += 1) {
        tot += 2.0 * path.getVint( path.bead(iPart,iBead) , path.bead(jPart,iBead) );
      }
    }
  }

  return tot;
}

// Get External Potential Energy
double Energy::getVEext()
{
  double tot = 0.0;

  if(path.trap) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
      for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {\
        tot += dot( path.bead(iPart,iBead) -> r , path.bead(iPart,iBead) -> r );
      }
    }
  }

  return path.halfOmega2 * tot;
}

// Get Nodal Energy
double Energy::getNE()
{
  if(!path.useNodeDist) return 0;

  for (unsigned int iPart = 0; iPart < path.nPart; iPart++) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead++) {
      path.updateNodeDistance(iPart,iBead);
    }
  }

  double xi, tot;
  tot = 0.0;
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
      double nD1 = path.bead(iPart,iBead) -> nDist;
      double nD2 = path.bead(iPart,iBead) -> next -> nDist;
      double nD1nD2 = nD1 * nD2;
      if (!nD1) {
        nD1nD2 = nD2 * nD2;
      } else if (!nD2) {
        nD1nD2 = nD1 * nD1;
      }
      xi = nD1nD2*path.oneOverLamTau;
      tot += xi/(path.tau*expm1(xi));
    }
  }

  return tot;

}
