#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder");
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  maxLevel = in.getAttribute<int>("maxLevel");
  int offset = 0;
  for (unsigned int iS=0; iS<path.nSpecies; iS+=1) {
    if (path.speciesList[iS]->name == speciesA) {
      iSpeciesA = iS;
      offsetA = offset;
    } else if (path.speciesList[iS]->name == speciesB) {
      iSpeciesB = iS;
      offsetB = offset;
    }
    offset += path.speciesList[iS]->nPart;
  }

  out.Write("/Actions/"+name+"/nImages", nImages);
  out.Write("/Actions/"+name+"/nOrder", nOrder);
  out.Write("/Actions/"+name+"/iSpeciesA", iSpeciesA);
  out.Write("/Actions/"+name+"/iSpeciesB", iSpeciesB);
}

RealType PairAction::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
      for (int iB=0; iB<path.nBead; iB+=1) {
        path.Dr(path(iP,iB),path(jP,iB),dr);
        RealType sum = 0.;
        for (int iD=0; iD<path.nD; iD++) {
          for (int image=-nImages; image<=nImages; image++) {
            RealType dist = dr(iD) + (RealType)image*path.L;
            sum += -0.5/dist;
          }
        }
        tot -= sum;
      }
    }
  }

  return tot;
}

RealType PairAction::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  if (level > maxLevel)
    return 0.;
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  Tvector dr;
  for (int iP=0; iP<particles.size(); ++iP) {
    for (int iB=b0; iB<b1; iB+=skip) {
      int iS = 0;
      if (path(iP,iB)->species.name == speciesA)
        iS = iSpeciesB;
      else if (path(iP,iB)->species.name == speciesB)
        iS = iSpeciesA;
      else
        cerr << "ERROR: Unrecognized species in PairAction action." << endl;
      for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
        path.Dr(path(iP,iB),path(jP,iB),dr);
        RealType gaussProd = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          RealType gaussSum = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            RealType dist = dr(iD) + (RealType)image*path.L;
            gaussSum += exp(-levelTau/dist);
          }
          gaussProd *= gaussSum;
        }
        tot -= log(gaussProd);
      }
    }
  }

  return tot;
}

void PairAction::Write()
{

}
