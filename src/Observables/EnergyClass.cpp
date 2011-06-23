#include "EnergyClass.h"

void Energy::Accumulate( const int pType )
{
  KE = path.getKE(); // Kinetic Energy
  VE = path.getVE(); // Potential Energy  
  NE = path.getNE(); // Nodal Energy      
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
