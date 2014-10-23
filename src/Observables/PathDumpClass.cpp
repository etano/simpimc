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

    // Get positions
    Tmatrix pathPositions(path.nPart*path.nBead,path.nD);
    for (int iP=0; iP<path.nPart; ++iP)
      for (int iB=0; iB<path.nBead; ++iB)
        for (int iD=0; iD<path.nD; ++iD)
          pathPositions(iP*path.nBead + iB,iD) = path(iP,iB)->r(iD);

    // Get permutation
    Tmatrix pathPermutation(path.nPart,2);
    for (int iP=0; iP<path.nPart; ++iP) {
      pathPermutation(iP,0) = path(iP,0)->prev->p;
      pathPermutation(iP,1) = path(iP,path.nBead-1)->next->p;
    }

    // Write to file
    if (firstTime) {
      firstTime = 0;
      out.Write("/Observables/"+name+"/nDump", nDump);
      out.CreateExtendableDataSet("/Observables/"+name+"/", "positions", pathPositions);
      out.CreateExtendableDataSet("/Observables/"+name+"/", "permutation", pathPermutation);
    } else {
      out.Rewrite("/Observables/"+name+"/nDump", nDump);
      out.AppendDataSet("/Observables/"+name+"/", "positions", pathPositions);
      out.AppendDataSet("/Observables/"+name+"/", "permutation", pathPermutation);
    }

    Reset();
  }

  nWriteCalls += 1;
}
