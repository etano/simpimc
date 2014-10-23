#include "PairCorrelationClass.h"

void PairCorrelation::Init(Input &in)
{
  // Read in species info
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  path.GetSpeciesInfo(speciesA, iSpeciesA, offsetA);
  path.GetSpeciesInfo(speciesB, iSpeciesB, offsetB);

  // Read in grid info
  RealType rMin = in.getAttribute<RealType>("rMin",0.);
  RealType rMax = in.getAttribute<RealType>("rMax",path.L/2.);
  int nR = in.getAttribute<RealType>("nR",100);
  gr.x.CreateGrid(rMin,rMax,nR);
  gr.y.zeros(nR);

  // Compute rs
  Tvector rs(gr.x.nR-1);
  for (int i=0; i<gr.x.nR-1; i++) {
    RealType r1 = gr.x(i);
    RealType r2 = gr.x(i+1);
    if (path.nD == 3)
      rs(i) = 0.75 * (r2*r2*r2*r2-r1*r1*r1*r1)/(r2*r2*r2-r1*r1*r1);
    else if (path.nD == 2)
      rs(i) = (r2*r2*r2-r1*r1*r1)/(r2*r2-r1*r1); // fixme: Not sure if 2D and 1D are correct here
    else if (path.nD == 1)
      rs(i) = 0.5*(r2-r1);
  }

  // Write things to file
  out.Write("Observables/"+name+"/speciesA", speciesA);
  out.Write("Observables/"+name+"/speciesB", speciesB);
  out.Write("Observables/"+name+"/rMin", rMin);
  out.Write("Observables/"+name+"/rMax", rMax);
  out.Write("Observables/"+name+"/nR", nR);
  out.Write("/Observables/"+name+"/r", rs);

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
  Tvector dr(path.nD);
  // Homogeneous
  if (iSpeciesA == iSpeciesB) {
    for (int iB=0; iB<path.nBead; ++iB) {
      for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart-1; ++iP) {
        for (int jP=iP+1; jP<offsetA+path.speciesList[iSpeciesA]->nPart; ++jP) {
          path.Dr(path(iP,iB),path(jP,iB),dr);
          int i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.nR)
            gr.y(i) = gr.y(i) + 1.;
        }
      }
    }
  // Homologous
  } else {
    for (int iB=0; iB<path.nBead; ++iB) {
      for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
          path.Dr(path(iP,iB),path(jP,iB),dr);
          int i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.nR)
            gr.y(i) = gr.y(i) + 1.;
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
    RealType vol = path.vol;
    RealType norm;
    if (iSpeciesA == iSpeciesB)
      norm = 0.5*nMeasure*NA*(NA-1.)*path.nBead/vol;
    else
      norm = 0.5*nMeasure*NA*NB*path.nBead/vol;
    for (int i=0; i<gr.x.nR; i++) {
      RealType r1 = gr.x(i);
      RealType r2 = (i<(gr.x.nR-1)) ? gr.x(i+1):(2.0*gr.x(i)-gr.x(i-1));
      RealType r = 0.5*(r1+r2);
      RealType binVol;
      if (path.nD == 3)
        binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      else if (path.nD == 2)
        binVol = M_PI * (r2*r2-r1*r1);
      else if (path.nD == 1)
        binVol = r2-r1;
      gr.y(i) = gr.y(i)/(binVol*norm);
    }

    // Write to file
    if (firstTime) {
      firstTime = 0;

      out.CreateExtendableDataSet("/Observables/"+name+"/", "gr", gr.y);
    } else {
      out.AppendDataSet("/Observables/"+name+"/", "gr", gr.y);
    }

    Reset();
  }
}
