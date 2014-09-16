#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder");
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  maxLevel = in.getAttribute<int>("maxLevel");
  useLongRange = in.getAttribute<int>("useLongRange",0);
  if (useLongRange) {
    kCut = in.getAttribute<RealType>("kCut");
    path.SetupKs(kCut);
  }
  GetOffset(speciesA,iSpeciesA,offsetA);
  GetOffset(speciesB,iSpeciesB,offsetB);

  string fileName = in.getAttribute<string>("file");
  ReadFile(fileName);

  //out.Write("Actions/"+name+"/file", fileName);
  //out.Write("Actions/"+name+"/nImages", nImages);
  //out.Write("Actions/"+name+"/nOrder", nOrder);
  //out.Write("Actions/"+name+"/speciesA", speciesA);
  //out.Write("Actions/"+name+"/speciesB", speciesB);
}

RealType PairAction::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
      if (iP != jP) {
        for (int iB=0; iB<path.nBead; iB+=1) {
          int jB = iB + 1;
          Tvector r, rP;
          path.Dr(path(iP,iB),path(jP,iB),r);
          path.Dr(path(iP,jB),path(jP,jB),rP);
          tot += CalcdUdBeta(r,rP,0);
        }
      }
    }
  }

  if (useLongRange) {
    tot += CalcdUdBetaLong();
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
      int jB = iB + skip;
      int iS = 0;
      int offset = 0;
      if (path(iP,iB)->species.name == speciesA) {
        iS = iSpeciesB;
        offset = offsetB;
      } else if (path(iP,iB)->species.name == speciesB) {
        iS = iSpeciesA;
        offset = offsetA;
      } else
        return 0.;
      for (int jP=offset; jP<offset+path.speciesList[iS]->nPart; ++jP) {
        if (iP != jP) {
          Tvector r, rP;
          path.Dr(path(iP,iB),path(jP,iB),r);
          path.Dr(path(iP,jB),path(jP,jB),rP);
          tot += CalcU(r,rP,level);
        }
      }
    }
  }

  if (useLongRange)
    tot += CalcULong(b0, b1, particles, level);

  return tot;
}

void PairAction::Write() {}

void PairAction::GetLimits(RealType &rMin, RealType &rMax, RealType &r, RealType &rP, NUgrid *g)
{
  rMax = g->end;
  if (r > rMax)
    r = rMax;
  if (rP > rMax)
    rP = rMax;
  rMin = g->start;
  if (rP < rMin)
    rP = rMin;
  if(r < rMin)
    r = rMin;
}

void PairAction::GetConstants(Tvector &rVec, Tvector &rPVec, RealType &r, RealType &rP, RealType &q, RealType &s, RealType &z, RealType &x, RealType &y)
{
  r = mag(rVec);
  rP = mag(rPVec);
  q = (r + rP)/2.;
  Tvector dr;
  path.Dr(rVec,rPVec,dr);
  s = mag(dr);
  z = r - rP;
  x = q + 0.5*s;
  y = q - 0.5*s;
}

