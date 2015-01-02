#include "KineticClass.h"

void Kinetic::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  species = in.getAttribute<string>("species");
  speciesList.push_back(species);
  cout << "Setting up kinetic action for " << species << "..." << endl;
  path.GetSpeciesInfo(species,iSpecies);
  nPart = path.speciesList[iSpecies]->nPart;
  i4LambdaTau = 1./(4.*path.speciesList[iSpecies]->lambda*path.tau);

  // Write things to file
  out.Write("Actions/"+name+"/nImages", nImages);
  out.Write("Actions/"+name+"/species", species);
}

double Kinetic::DActionDBeta()
{
  double tot = nPart*path.nBead*path.nD/(2.*path.tau); // Constant term
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iP=0; iP<nPart; iP++) {
    for (int iB=0; iB<path.nBead; iB++) {
      vec<double> numSum, gaussSum, dr(path.nD);
      path.Dr(path(iSpecies,iP,iB),path.GetNextBead(path(iSpecies,iP,iB),1),dr);
      double gaussProd = 1.;
      numSum.zeros(path.nD);
      gaussSum.zeros(path.nD);
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
  
  
  return tot;
}

double Kinetic::GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  double tot = 0.;
  vec<double> dr(path.nD);
  std::shared_ptr<Bead> beadA, beadB, beadC, beadF;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      double gaussProd, gaussSum, dist;
      beadA = path(iSpecies,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        beadB = path.GetNextBead(beadA,skip);
        path.Dr(beadA,beadB,dr);
        gaussProd = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          gaussSum = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            dist = dr(iD) + image*path.L;
            gaussSum += exp(-dist*dist*i4LambdaLevelTau);
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
