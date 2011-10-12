#include "RClass.h"

void R::Accumulate( const int pType )
{
  Rtemp.zeros();
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      Rtemp(iPart) += (path.bead(iPart,iBead) -> r)(0); // Add up R
    }
    Rtot(pType,iPart) += Rtemp(iPart); // Sum total R
  }  
}

void R::Output()
{
  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
      trace << Rtot(iType,iPart) * oneOverNbeadBlock << " ";
    }
  }
  trace << endl;

  Rtot.zeros();
}

void R::Print()
{
  std::cout << "\nR: ";
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    std::cout << endl << iPart << " : " << Rtemp(iPart); 
  }
  std::cout << endl;
}

void R::Stats()
{
  statsR( outputFile.c_str() , path.nType , nBlock , path.nPart , path.nD );
}
