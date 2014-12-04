#include "KineticClass.h"

void Kinetic::Init(Input &in)
{
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;

  // todo: spline maybe?

  normalization = path.nBead*path.nD/(2.*path.tau);
  out.Write("/Actions/"+name+"/nImages", nImages);
}

double Kinetic::DActionDBeta()
{
  double tot = 0.;
  for (int iS=0; iS<path.nSpecies; iS++) {
    double lambda = path.speciesList[iS]->lambda;
    if (!fequal(lambda,0.,1e-10)) {
      double i4LambdaTau = 1./(4.*lambda*path.tau);
      tot += path.speciesList[iS]->nPart*path.nBead*path.nD/(2.*path.tau); // Constant term
      #pragma omp parallel for collapse(2) reduction(+:tot)
      for (int iP=0; iP<path.speciesList[iS]->nPart; iP++) {
        for (int iB=0; iB<path.nBead; iB++) {
          vec<double> gaussSum, numSum, dr(path.nD);
          path.Dr(path(iS,iP,iB),path.GetNextBead(path(iS,iP,iB),1),dr);
          double gaussProd = 1.;
          gaussSum.zeros(path.nD);
          numSum.zeros(path.nD);
          for (int iD=0; iD<path.nD; iD++) {
            for (int image=-nImages; image<=nImages; image++) {
              double dist = dr(iD) + image*path.L;
              double dist2i4LambdaTau = dist*dist*i4LambdaTau;
              double expPart = path.fexp(-dist2i4LambdaTau);
              gaussSum(iD) += expPart;
              numSum(iD) += -(dist2i4LambdaTau/path.tau)*expPart;
            }
            gaussProd *= gaussSum(iD);
          }
          double scalarNumSum = 0.;
          for (int iD=0; iD<path.nD; iD++) {
            double numProd = 1.;
            for (int jD=0; jD<path.nD; jD++) {
              if (iD != jD)
                numProd *= gaussSum(jD);
              else
                numProd *= numSum(jD);
            }
            scalarNumSum += numProd;
          }
          tot += scalarNumSum/gaussProd;
        }
      }
    }
  }
  return tot;
}

double Kinetic::GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double levelTau = skip*path.tau;
  double tot = 0.;
  vec<double> dr(path.nD);
  Bead *beadA, *beadB, *beadC, *beadF;
  for (int p=0; p<particles.size(); ++p) {
    int iS = particles[p].first;
    int iP = particles[p].second;
    double lambda = path.speciesList[iS]->lambda;
    if (!fequal(lambda,0.,1e-10)) {
      double i4LambdaTau = 1./(4.*lambda*levelTau);
      double gaussProd, gaussSum, dist;
      beadA = path(iS,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        beadB = path.GetNextBead(beadA,skip);
        path.Dr(beadA,beadB,dr);
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

        beadA = beadB;
      }
    }
  }

  return tot;
}

void Kinetic::Write()
{

}
