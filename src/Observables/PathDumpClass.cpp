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
      Tmatrix pathPositions(nPart*path.nBead,path.nD);
      for (int iP=0; iP<nPart; ++iP)
        for (int iB=0; iB<path.nBead; ++iB)
          for (int iD=0; iD<path.nD; ++iD)
            pathPositions(iP*path.nBead + iB,iD) = path(iS,iP,iB)->r(iD);

      // Get permutation
      Tmatrix pathPermutation(nPart,2);
      for (int iP=0; iP<nPart; ++iP) {
        pathPermutation(iP,0) = path(iS,iP,0)->prev->p;
        pathPermutation(iP,1) = path(iS,iP,path.nBead-1)->next->p;
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
