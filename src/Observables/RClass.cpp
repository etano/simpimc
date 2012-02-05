#include "RClass.h"

void R::Accumulate( const int pType )
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      Rblock(pType,iPart) += (path.bead(iPart,iBead) -> r)(0); // Add up R
    }
  }
}

void R::Output()
{
  Rblock *= oneOverNbeadBlock;
  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
      trace << Rblock(iType,iPart) << " ";
    }
  }
  trace << endl;

  for (unsigned int iPart = 0; iPart < path.nPart; iPart++) {
    Rtot(iPart) += sum(Rblock.col(iPart));
  }

  trace << endl;
  Rblock.zeros();
}

void R::Print()
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    std::cout << iPart << " : " << Rtot(iPart)/nBlock << endl;
  }
}

void R::Stats()
{
  statsR( outputFile.c_str() , path.nType , nBlock , path.nPart , path.nD );
}
