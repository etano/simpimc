#include "EnergyClass.h"

void Energy::Accumulate( const int iType )
{
  KE(iType) += getKE(); // Add up total Kinetic energy
  VE(iType) += getVE(); // Add up total Potential energy
  NE(iType) += getNE(); // Add up total Nodal energy
  E(iType) = KE(iType) + VE(iType) + NE(iType); // Add up total energy
}

void Energy::Output()
{
  E *= oneOverNbeadBlock;
  KE *= oneOverNbeadBlock;
  VE *= oneOverNbeadBlock;
  NE *= oneOverNbeadBlock;

  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    trace << E(iType) << " " << KE(iType) << " " << VE(iType) << " " << NE(iType) << " ";
  }
  trace << sum(E) << " " << sum(KE) << " " << sum(VE) << " " << sum(NE) << endl;

  Etot += sum(E);
  KEtot += sum(KE);
  VEtot += sum(VE);
  NEtot += sum(NE);

  E.zeros();
  KE.zeros();
  VE.zeros();
  NE.zeros();
}

void Energy::Print()
{
  std::cout << "E: " << Etot/nBlock << endl
            << "KE: " << KEtot/nBlock << endl
            << "VE: " << VEtot/nBlock << endl
            << "NE: " << NEtot/nBlock << endl;
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
      for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
        tot += dot( path.bead(iPart,iBead) -> r , path.bead(iPart,iBead) -> r );
      }
    }
  }

  return path.halfOmega2 * path.onePlus3Tau2Omega2Over12 * tot;
}

// Get Nodal Energy
double Energy::getNE()
{
  if(!path.useNodeDist) return 0;

  //path.mode = 1;
  //for (unsigned int iPart = 0; iPart < path.nPart; iPart++) {
  //  for (unsigned int iBead = 0; iBead < path.nBead; iBead++) {
  //    path.updateNodeDistance(path.bead(iPart,iBead));
  //  }
  //}

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
      //xi = nD1nD2*path.oneOverLamTau;
      xi = 0.5*nD1nD2*path.oneOverLamTau;
      tot += xi/(path.tau*expm1(xi));
      if (std::isnan(tot))
        std::cerr << tot << endl;
    }
  }

  return tot;

}
