#include "PairCorrelationClass.h"

void PairCorrelation::Init(Input &in)
{
  // Read in species info
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  path.GetSpeciesInfo(speciesA, iSpeciesA);
  path.GetSpeciesInfo(speciesB, iSpeciesB);

  // Read in grid info
  double rMin = in.getAttribute<double>("rMin",0.);
  double rMax = in.getAttribute<double>("rMax",path.L/2.);
  int nR = in.getAttribute<double>("nR",1000);
  gr.x.CreateGrid(rMin,rMax,nR);
  gr.y.zeros(nR);

  // Compute rs
  vec<double> rs(gr.x.nR-1);
  for (int i=0; i<gr.x.nR-1; i++) {
    double r1 = gr.x(i);
    double r2 = gr.x(i+1);
    if (path.nD == 3)
      rs(i) = 0.75 * (r2*r2*r2*r2-r1*r1*r1*r1)/(r2*r2*r2-r1*r1*r1);
    else if (path.nD == 2)
      rs(i) = (r2*r2*r2-r1*r1*r1)/(r2*r2-r1*r1); // fixme: Not sure if 2D and 1D are correct here
    else if (path.nD == 1)
      rs(i) = 0.5*(r2-r1);
  }

  // Write things to file
  out.Write(prefix+"/speciesA", speciesA);
  out.Write(prefix+"/speciesB", speciesB);
  out.Write(prefix+"/rMin", rMin);
  out.Write(prefix+"/rMax", rMax);
  out.Write(prefix+"/nR", nR);
  out.Write(prefix+"/x", rs);

  Reset();
}

void PairCorrelation::Reset()
{
  nMeasure = 0;
  gr.y.zeros();
}

void PairCorrelation::Accumulate()
{
  path.SetMode(1);
  vec<double> dr(path.nD);
  // Homogeneous
  if (iSpeciesA == iSpeciesB) {
    for (int iB=0; iB<path.nBead; ++iB) {
      for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
        for (int jP=iP+1; jP<path.speciesList[iSpeciesA]->nPart; ++jP) {
          path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesA,jP,iB),dr);
          int i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.nR)
            gr.y(i) = gr.y(i) + 1.*path.importance_weight;
        }
      }
    }
  // Homologous
  } else {
    for (int iB=0; iB<path.nBead; ++iB) {
      for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (int jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
          path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesB,jP,iB),dr);
          int i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.nR)
            gr.y(i) = gr.y(i) + 1.*path.importance_weight;
        }
      }
    }
  }

  nMeasure += 1;
}

void PairCorrelation::Write()
{
  if (nMeasure > 0) {
    // Normalize histograms
    int NA = path.speciesList[iSpeciesA]->nPart;
    int NB = path.speciesList[iSpeciesB]->nPart;
    double vol = path.vol;
    double norm;
    if (iSpeciesA == iSpeciesB)
      norm = 0.5*nMeasure*NA*(NA-1.)*path.nBead/vol;
    else
      norm = nMeasure*NA*NB*path.nBead/vol;
    for (int i=0; i<gr.x.nR; i++) {
      double r1 = gr.x(i);
      double r2 = (i<(gr.x.nR-1)) ? gr.x(i+1):(2.*gr.x(i)-gr.x(i-1));
      double r = 0.5*(r1+r2);
      double binVol;
      if (path.nD == 3)
        binVol = 4.*M_PI/3. * (r2*r2*r2-r1*r1*r1);
      else if (path.nD == 2)
        binVol = M_PI * (r2*r2-r1*r1);
      else if (path.nD == 1)
        binVol = r2-r1;
      //gr.y(i) = gr.y(i)/(binVol*norm);
      gr.y(i) = gr.y(i)/(norm);
    }

    // Write to file
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet("/"+prefix, "y", gr.y);
    } else {
      out.AppendDataSet("/"+prefix, "y", gr.y);
    }

    Reset();
  }
}
