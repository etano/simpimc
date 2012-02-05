#include "R2Class.h"

void R2::Accumulate( const int pType )
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1) {
      R2block(pType,iPart) += dot( path.bead(iPart,iBead)->r , path.bead(iPart,iBead)->r ); // Add up R2
    }
  }
}

void R2::Output()
{
  R2block *= oneOverNbeadBlock;
  for (unsigned int iType = 0; iType < path.nType; iType += 1) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
      trace << R2block(iType,iPart) << " ";
    }
  }

  for (unsigned int iPart = 0; iPart < path.nPart; iPart++) {
    R2tot(iPart) += sum(R2block.col(iPart));
  }

  trace << endl;
  R2block.zeros();
}

void R2::Print()
{
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1) {
    std::cout << iPart << " : " << R2tot(iPart)/nBlock << endl;
  }
}

void R2::Stats()
{
  statsRR( outputFile.c_str() , path.nType , nBlock , path.nPart );
}
