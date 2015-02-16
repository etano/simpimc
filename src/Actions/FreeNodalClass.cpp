#include "FreeNodalClass.h"

// Initialize parameters
void FreeNodal::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  species = in.getAttribute<string>("species");
  speciesList.push_back(species);
  cout << "Setting up nodal action for " << species << "..." << endl;
  maxLevel = in.getAttribute<int>("maxLevel",0);
  path.GetSpeciesInfo(species,iSpecies);
  nPart = path.speciesList[iSpecies]->nPart;
  i4LambdaTau = 1./(4.*path.speciesList[iSpecies]->lambda*path.tau);

  // Write things to file
  out.Write("Actions/"+name+"/nImages", nImages);
  out.Write("Actions/"+name+"/species", species);
  out.Write("Actions/"+name+"/maxLevel", maxLevel);

  // Setup splines
  SetupSpline();

  // Set up determinants
  rho_F.set_size(path.nBead);
  rho_F_c.set_size(path.nBead);

  // Test initial configuration
  std::vector< std::pair<int,int> > particles;
  for (int iP=0; iP<nPart; ++iP)
    particles.push_back(std::make_pair(iSpecies,iP));
  int initGood = 1;
  path.SetMode(1);
  if (GetAction(0, path.nBead, particles, 0) == 1.e100) {
    cout << "Warning: initializing with broken nodes!" << endl;
    initGood = 0;
  }
  out.Write("Actions/"+name+"/initGood", initGood);
}

// Create a spline for each possible sliceDiff
void FreeNodal::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  r_grid.start = -path.L/2.;
  r_grid.end = path.L/2.;
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  int nSpline = path.nBead/2 + (path.nBead%2) + 1;
  rho_free_r_splines.set_size(nSpline);

  // Create splines
  for (int iSpline=0; iSpline<nSpline; ++iSpline) {
    vec<double> rho_free_r(r_grid.num);
    double t_i4LambdaTau = i4LambdaTau/(iSpline+1);

    // Make rho_0
    double rho0 = 0.;
    for (int image=-nImages; image<=nImages; image++)
      rho0 += path.fexp(-path.L*path.L*t_i4LambdaTau);
    double logRho0 = -log(rho0);

    // Make rho_free
    for (int i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      rho_free_r(i) = 0.;
      for (int image=-nImages; image<=nImages; image++) {
        if (image != 0) {
          double t_r = r + image*path.L;
          rho_free_r(i) += path.fexp(-(t_r*t_r - r*r)*t_i4LambdaTau);
        }
      }
      rho_free_r(i) = log1p(min(10.,rho_free_r(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines(iSpline) = rho_free_r_spline;
  }
}

// Evaluate \rho_{ij}
double FreeNodal::GetGij(vec<double>& r, int sliceDiff)
{
  double gaussProd = 1.;
  for (int iD=0; iD<path.nD; iD++) {
    double gaussSum;
    eval_UBspline_1d_d(rho_free_r_splines(sliceDiff-1),r(iD),&gaussSum);
    gaussSum = exp(0.9999*gaussSum);
    gaussSum *= exp(-(r(iD)*r(iD)*i4LambdaTau/sliceDiff));
    gaussProd *= gaussSum;
  }
  return gaussProd;
}
