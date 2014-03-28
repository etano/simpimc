#include "CoulombClass.h"

void Coulomb::Init(Input &in)
{
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder");
  Z1Z2 = in.getAttribute<RealType>("Z1Z2");
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  maxLevel = in.getAttribute<int>("maxLevel");
  GetOffset(speciesA,iSpeciesA,offsetA);
  GetOffset(speciesB,iSpeciesB,offsetB);

  string fileName = "hello";

  out.Write("/Actions/"+name+"/file", fileName);
  out.Write("/Actions/"+name+"/nImages", nImages);
  out.Write("/Actions/"+name+"/nOrder", nOrder);
  out.Write("/Actions/"+name+"/Z1Z2", Z1Z2);
  out.Write("/Actions/"+name+"/speciesA", speciesA);
  out.Write("/Actions/"+name+"/speciesB", speciesB);
}

RealType Coulomb::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
      for (int iB=0; iB<path.nBead; iB+=1) {
        path.Dr(path(iP,iB),path(jP,iB),dr);
        tot += Z1Z2/sqrt(dot(dr,dr));
        //RealType sum = 0.;
        //for (int iD=0; iD<path.nD; iD++) {
        //  for (int image=-nImages; image<=nImages; image++) {
        //    RealType dist = dr(iD) + (RealType)image*path.L;
        //    sum += -0.5*Z1Z2/dist;
        //  }
        //}
        //tot -= sum;
      }
    }
  }

  return tot;
}

RealType Coulomb::GetAction(int b0, int b1, vector<int> &particles, int level)
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
      int offset = 0;
      if (path(iP,iB)->species.name == speciesA) {
        iS = iSpeciesB;
        offset = offsetB;
      } else if (path(iP,iB)->species.name == speciesB) {
        iS = iSpeciesA;
        offset = offsetA;
      } else
        cerr << "ERROR: Unrecognized species in Coulomb action." << endl;
      for (int jP=offset; jP<offset+path.speciesList[iS]->nPart; ++jP) {
        path.Dr(path(iP,iB),path(jP,iB),dr);
        tot += levelTau*Z1Z2/sqrt(dot(dr,dr));
        //RealType gaussProd = 1.;
        //for (int iD=0; iD<path.nD; iD++) {
        //  RealType gaussSum = 0.;
        //  for (int image=-nImages; image<=nImages; image++) {
        //    RealType dist = dr(iD) + (RealType)image*path.L;
        //    gaussSum += exp(-levelTau*Z1Z2/dist);
        //  }
        //  gaussProd *= gaussSum;
        //}
        //tot -= log(gaussProd);
      }
    }
  }

  return tot;
}

void Coulomb::Write()
{

}
