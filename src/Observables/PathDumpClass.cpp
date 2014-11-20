#include "PathDumpClass.h"

void PathDump::Init(Input &in)
{
  nDump = 0;

  Reset();
}

void PathDump::Reset()
{
  nWriteCalls = 0;
}

void PathDump::Accumulate() {}

void PathDump::Write()
{
  if (!(nWriteCalls % skip)) {
    nDump += 1;

    // Loop through species
    for (int iS=0; iS<path.nSpecies; ++iS) {
       // Get species info
      path.GetSpeciesInfo(path.speciesList[iS]->name,iS);
      int nPart = path.speciesList[iS]->nPart;

      // Get positions
      Tcube pathPositions(path.nD,path.nBead,nPart);
      for (int iP=0; iP<nPart; ++iP)
        for (int iB=0; iB<path.nBead; ++iB)
          for (int iD=0; iD<path.nD; ++iD)
            pathPositions(iD,iB,iP) = path(iS,iP,iB)->r(iD);

      // Get permutation
      Tmatrix pathPermutation(2,nPart);
      for (int iP=0; iP<nPart; ++iP) {
        pathPermutation(0,iP) = path(iS,iP,0)->prev->p;
        pathPermutation(1,iP) = path(iS,iP,path.nBead-1)->next->p;
      }

      // Write to file
      if (firstTime) {
        out.CreateGroup(prefix+path.speciesList[iS]->name+"/");
        out.Write(prefix+path.speciesList[iS]->name+"/nDump", nDump);
        out.CreateExtendableDataSet(prefix+path.speciesList[iS]->name+"/", "positions", pathPositions);
        out.CreateExtendableDataSet(prefix+path.speciesList[iS]->name+"/", "permutation", pathPermutation);
      } else {
        out.Rewrite(prefix+path.speciesList[iS]->name+"/nDump", nDump);
        out.AppendDataSet(prefix+path.speciesList[iS]->name+"/", "positions", pathPositions);
        out.AppendDataSet(prefix+path.speciesList[iS]->name+"/", "permutation", pathPermutation);
      }

    }

    if (firstTime)
      firstTime = 0;

    Reset();
  }

  nWriteCalls += 1;
}
