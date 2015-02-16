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

  // Setup spline
  SetupSpline();
}

// Create a spline for each possible sliceDiff
// TODO: Combine with free nodal action splines
void Kinetic::SetupSpline()
{
  // Setup grid
  Ugrid r2_grid;
  r2_grid.start = 1.e-8;
  r2_grid.end = path.L/2.;
  r2_grid.num = 10000;
  double dr = (r2_grid.end - r2_grid.start)/(r2_grid.num - 1);

  // Resize spline field
  int nSpline = path.nBead/2 + (path.nBead%2) + 1;
  rho_free_r2_splines.set_size(nSpline);

  // Create splines
  for (int iSpline=0; iSpline<nSpline; ++iSpline) {
    vec<double> rho_free_r2(r2_grid.num), num_sum_r2(r2_grid.num);
    double t_i4LambdaTau = i4LambdaTau/(iSpline+1);

    // Make rho_free
    for (int i=0; i<r2_grid.num; ++i) {
      double r2 = r2_grid.start + i*dr;
      rho_free_r2(i) = 0.;
      if (iSpline == 0)
        num_sum_r2(i) = 0.;
      double r2i4LambdaTau = r2*t_i4LambdaTau;
      double L2 = path.L*path.L;
      for (int image=-nImages; image<=nImages; image++) {
        double t_r2 = r2 + 2*sqrt(r2)*image*path.L + image*image*L2;
        double t_r2i4LambdaTau = t_r2*t_i4LambdaTau;
        if (image != 0) {
          double expPart = path.fexp(r2i4LambdaTau-t_r2i4LambdaTau);
          rho_free_r2(i) += expPart;
          if (iSpline == 0)
            num_sum_r2(i) += (t_r2/r2)*expPart;
        }
      }
      rho_free_r2(i) = log1p(min(10.,rho_free_r2(i)));
      num_sum_r2(i) = log1p(min(10.,num_sum_r2(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r2_spline = create_UBspline_1d_d(r2_grid, xBC, rho_free_r2.memptr());
    rho_free_r2_splines(iSpline) = rho_free_r2_spline;
    if (iSpline == 0)
      num_sum_r2_spline = create_UBspline_1d_d(r2_grid, xBC, num_sum_r2.memptr());
  }
}

double Kinetic::GetGaussSum(const double &r, const int sliceDiff)
{
  double gaussSum;
  double r2 = r*r;
  eval_UBspline_1d_d(rho_free_r2_splines(sliceDiff-1),r2,&gaussSum);
  gaussSum = exp(0.9999*gaussSum);
  gaussSum *= exp(-(r2*i4LambdaTau/sliceDiff));
  return gaussSum;
}

double Kinetic::GetNumSum(const double &r)
{
  double numSum;
  double r2 = r*r;
  eval_UBspline_1d_d(num_sum_r2_spline,r2,&numSum);
  numSum = exp(0.9999*numSum);
  numSum *= -(r2*i4LambdaTau/path.tau)*exp(-(r2*i4LambdaTau));
  return numSum;
}

double Kinetic::DActionDBeta()
{
  double tot = nPart*path.nBead*path.nD/(2.*path.tau); // Constant term
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iP=0; iP<nPart; iP++) {
    for (int iB=0; iB<path.nBead; iB++) {
      vec<double> numSum, gaussSum, dr(path.nD);
      path.Dr(path(iSpecies,iP,iB),path.GetNextBead(path(iSpecies,iP,iB),1),dr);
      numSum.zeros(path.nD);
      gaussSum.zeros(path.nD);
      double gaussProd = 1.;
      for (int iD=0; iD<path.nD; iD++) {
        numSum(iD) = GetNumSum(dr(iD));
        gaussSum(iD) = GetGaussSum(dr(iD),1);
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
  std::shared_ptr<Bead> beadA, beadB, beadF;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      beadA = path(iSpecies,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        beadB = path.GetNextBead(beadA,skip);
        path.Dr(beadA,beadB,dr);
        double gaussProd = 1;
        for (int iD=0; iD<path.nD; iD++)
          gaussProd *= GetGaussSum(dr(iD),skip);
        tot -= log(gaussProd);

        beadA = beadB;
      }
    }
  }

  return tot;
}

vec<double> Kinetic::GetActionGradient(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  vec<double> tot;
  tot.zeros(path.nD);
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
        beadB = path.GetPrevBead(beadA,skip);
        path.Dr(beadB,beadA,dr);
        tot -= dr;
        beadC = path.GetNextBead(beadA,skip);
        path.Dr(beadA,beadC,dr);
        tot += dr;
        beadA = beadC;
      }
    }
  }

  return 2.*i4LambdaLevelTau*tot;
}

double Kinetic::GetActionLaplacian(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  double tot = 0.;
  vec<double> dr(path.nD);
  std::shared_ptr<Bead> beadA, beadF;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      double gaussProd, gaussSum, dist;
      beadA = path(iSpecies,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        tot += path.nD*4.*i4LambdaLevelTau;
        beadA = path.GetNextBead(beadA,skip);
      }
    }
  }

  return tot;
}

void Kinetic::Write()
{

}
