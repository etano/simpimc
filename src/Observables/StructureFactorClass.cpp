#include "StructureFactorClass.h"

void StructureFactor::Init(Input &in)
{
  // Read in species info
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  path.GetSpeciesInfo(speciesA, iSpeciesA, offsetA);
  path.GetSpeciesInfo(speciesB, iSpeciesB, offsetB);
  kCut = in.getAttribute<RealType>("kCut", path.kC);

  // Resize
  path.SetupKs(kCut);
  sk.zeros(path.ks.size());

  // Write things to file
  out.Write("Observables/"+name+"/speciesA", speciesA);
  out.Write("Observables/"+name+"/speciesB", speciesB);
  out.Write("Observables/"+name+"/kCut", kCut);
  out.CreateExtendableDataSet("/Observables/"+name+"/", "ks", path.magKs[0]);
  for (int iK=1; iK<path.magKs.size(); ++iK)
    out.AppendDataSet("/Observables/"+name+"/", "ks", path.magKs[iK]);

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
  for (int iK=0; iK<path.kIndices.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (int iB=0; iB<path.nBead; ++iB) {
        sk(iK) += cmag2(path.rhoK(iB,iSpeciesA)(iK),path.rhoK(iB,iSpeciesB)(iK));
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
    int NA = path.speciesList[iSpeciesA]->nPart;
    int NB = path.speciesList[iSpeciesB]->nPart;
    RealType norm = nMeasure*path.nBead*NA*NB;
    sk = sk/norm;

    // Write to file
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet("/Observables/"+name+"/", "sk", sk);
    } else {
      out.AppendDataSet("/Observables/"+name+"/", "sk", sk);
    }

    Reset();
  }
}
