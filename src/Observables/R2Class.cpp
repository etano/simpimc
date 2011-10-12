#include "R2Class.h"

void R2::Accumulate( const int pType )
{
  R2temp.zeros();
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      R2temp(iPart) += dot( path.bead(iPart,iBead) -> r , path.bead(iPart,iBead) -> r ); // Add up R2
    }
    R2tot(pType,iPart) += R2temp(iPart); // Sum total R2
  }  
}

void R2::Output()
{
  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
      trace << R2tot(iType,iPart) * oneOverNbeadBlock << " ";
    }
  }
  trace << endl;

  R2tot.zeros();
}

void R2::Print()
{
  std::cout << "\nR2: ";
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    std::cout << endl << iPart << " : " << R2temp(iPart); 
  }
  std::cout << endl;
}

void R2::Stats()
{
  statsRR( outputFile.c_str() , path.nType , nBlock , path.nPart );
}
