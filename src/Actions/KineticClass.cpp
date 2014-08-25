#include "KineticClass.h"

void Kinetic::Init(Input &in)
{
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;
  //out.Write("/Actions/"+name+"/nImages", nImages);
}

RealType Kinetic::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=0; iP<path.nPart; iP+=1) {
    RealType lambda = path.speciesList[iP]->lambda;
    if (lambda != 0) {
      RealType i4LambdaTau = 1./(4.*lambda*path.tau);
      for (int iB=0; iB<path.nBead; iB+=1) {
        path.Dr(path(iP,iB),path(iP,iB)->next,dr);
        RealType gaussProd = 1.;
        Tvector gaussSum, numSum;
        gaussSum.zeros(path.nD);
        numSum.zeros(path.nD);
        for (int iD=0; iD<path.nD; iD++) {
          for (int image=-nImages; image<=nImages; image++) {
            RealType dist = dr(iD) + (RealType)image*path.L;
            RealType dist2i4LambdaTau = dist*dist*i4LambdaTau;
            RealType expPart = exp(-dist2i4LambdaTau);
            gaussSum(iD) += expPart;
            numSum(iD) += -(dist2i4LambdaTau/path.tau)*expPart;
          }
          gaussProd *= gaussSum(iD);
        }
        RealType scalarNumSum = 0.;
        for (int iD=0; iD<path.nD; iD++) {
          Tvector numProd;
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
  return tot;
}

RealType Kinetic::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  Tvector dr;
  for (int iP=0; iP<particles.size(); ++iP) {
    RealType lambda = path.speciesList[particles[iP]]->lambda;
    if (lambda != 0) {
      RealType i4LambdaTau = 1./(4.*lambda*levelTau);
      for (int iB=b0; iB<b1; iB+=skip) {
        path.Dr(path(iP,iB),path(iP,iB)->nextB(skip),dr);
        RealType gaussProd = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          RealType gaussSum = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            RealType dist = dr(iD) + (RealType)image*path.L;
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
