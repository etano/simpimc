#include "KineticClass.h"

void Kinetic::Init(Input &in)
{
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;

  normalization = path.nBead*path.nD/(2.*path.tau);
  out.Write("/Actions/"+name+"/nImages", nImages);
}

RealType Kinetic::DActionDBeta()
{
  RealType tot = 0.;
  Tvector gaussSum(path.nD), numSum(path.nD), numProd(path.nD), dr(path.nD);
  for (int iS=0; iS<path.nSpecies; iS++) {
    RealType lambda = path.speciesList[iS]->lambda;
    if (!fequal(lambda,0.,1e-10)) {
      RealType i4LambdaTau = 1./(4.*lambda*path.tau);
      RealType gaussProd, dist, dist2i4LambdaTau, expPart, scalarNumSum;
      int offset;
      GetOffset(path.speciesList[iS]->name,iS,offset);
      for (int iP=offset; iP<offset+path.speciesList[iS]->nPart; iP++) {
        for (int iB=0; iB<path.nBead; iB+=1) {
          path.Dr(path(iP,iB),path(iP,iB)->next,dr);
          gaussProd = 1.;
          gaussSum.zeros();
          numSum.zeros();
          for (int iD=0; iD<path.nD; iD++) {
            for (int image=-nImages; image<=nImages; image++) {
              dist = dr(iD) + image*path.L;
              dist2i4LambdaTau = dist*dist*i4LambdaTau;
              expPart = exp(-dist2i4LambdaTau);
              gaussSum(iD) += expPart;
              numSum(iD) += -(dist2i4LambdaTau/path.tau)*expPart;
            }
            gaussProd *= gaussSum(iD);
          }
          scalarNumSum = 0.;
          for (int iD=0; iD<path.nD; iD++) {
            numProd.ones(path.nD);
            for (int jD=0; jD<path.nD; jD++) {
              if (iD != jD)
                numProd(iD) *= gaussSum(jD);
              else
                numProd(iD) *= numSum(jD);
            }
            scalarNumSum += numProd(iD);
          }
          tot += scalarNumSum/gaussProd;
        }
        tot += path.nBead*path.nD/(2.*path.tau);
      }
    }
  }
  return tot;
}

RealType Kinetic::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  Tvector dr(path.nD);
  for (int p=0; p<particles.size(); ++p) {
    int iP = particles[p];
    int iS = path(iP,b0)->s;
    RealType lambda = path.speciesList[iS]->lambda;
    if (!fequal(lambda,0.,1e-10)) {
      RealType i4LambdaTau = 1./(4.*lambda*levelTau);
      RealType gaussProd, gaussSum, dist;
      for (int iB=b0; iB<b1; iB+=skip) {
        path.Dr(path(iP,iB),path(iP,iB)->nextB(skip),dr);
        gaussProd = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          gaussSum = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            dist = dr(iD) + image*path.L;
            gaussSum += exp(-dist*dist*i4LambdaTau);
          }
          gaussProd *= gaussSum;
        }
        tot -= log(gaussProd);
      }
    }
  }

  return tot;
}

void Kinetic::Write()
{

}
