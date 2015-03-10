#include "StructureFactorClass.h"

void StructureFactor::Init(Input &in)
{
  // Read in species info
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  path.GetSpeciesInfo(speciesA, iSpeciesA);
  path.GetSpeciesInfo(speciesB, iSpeciesB);
  kCut = in.getAttribute<double>("kCut", path.kC);

  // Resize
  path.SetupKs(kCut);
  sk.zeros(path.ks.size());

  // Write things to file
  out.Write(prefix+"/speciesA", speciesA);
  out.Write(prefix+"/speciesB", speciesB);
  out.Write(prefix+"/kCut", kCut);
  out.CreateExtendableDataSet("/"+prefix+"/", "x", path.magKs[0]);
  for (uint iK=1; iK<path.magKs.size(); ++iK)
    out.AppendDataSet("/"+prefix+"/", "x", path.magKs[iK]);

  Reset();
}

void StructureFactor::Reset()
{
  nMeasure = 0;
  sk.zeros();
}

void StructureFactor::Accumulate()
{
  path.SetMode(1);
  for (uint iK=0; iK<path.kIndices.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (uint iB=0; iB<path.nBead; ++iB) {
        sk(iK) += path.sign*path.importance_weight*cmag2(path.rhoK(iB,iSpeciesA)(iK),path.rhoK(iB,iSpeciesB)(iK));
      }
    }
  }

  //if (iSpeciesB != iSpeciesA)
  //  sk = 2.*sk;

  nMeasure += 1;
}

void StructureFactor::Write()
{
  if (nMeasure > 0) {
    // Normalize histograms
    uint NA = path.speciesList[iSpeciesA]->nPart;
    uint NB = path.speciesList[iSpeciesB]->nPart;
    double norm = nMeasure*path.nBead*NA*NB;
    sk = sk/norm;

    // Write to file
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet("/"+prefix+"/", "y", sk);
    } else {
      out.AppendDataSet("/"+prefix+"/", "y", sk);
    }

    Reset();
  }
}
