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
    vec<double> rho_free_r2(r2_grid.num);
    double t_i4LambdaTau = i4LambdaTau/(iSpline+1);

    // Make rho_free
    for (int i=0; i<r2_grid.num; ++i) {
      double r2 = r2_grid.start + i*dr;
      rho_free_r2(i) = 0.;
      double r2i4LambdaTau = r2*t_i4LambdaTau;
      double L2 = path.L*path.L;
      for (int image=-nImages; image<=nImages; image++) {
        double t_r2 = r2 + 2*sqrt(r2)*image*path.L + image*image*L2;
        double t_r2i4LambdaTau = t_r2*t_i4LambdaTau;
        if (image != 0)
          rho_free_r2(i) += path.fexp(r2i4LambdaTau-t_r2i4LambdaTau);
      }
      rho_free_r2(i) = log1p(min(10.,rho_free_r2(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r2_spline = create_UBspline_1d_d(r2_grid, xBC, rho_free_r2.memptr());
    rho_free_r2_splines(iSpline) = rho_free_r2_spline;
  }
}

// Evaluate \rho_{ij}
double FreeNodal::GetGij(vec<double>& r, int sliceDiff)
{
  double gaussProd = 1.;
  for (int iD=0; iD<path.nD; iD++) {
    double gaussSum;
    double r2 = r(iD)*r(iD);
    eval_UBspline_1d_d(rho_free_r2_splines(sliceDiff-1),r2,&gaussSum);
    gaussSum = exp(0.9999*gaussSum);
    gaussSum *= exp(-(r2*i4LambdaTau/sliceDiff));
    gaussProd *= gaussSum;
  }
  return gaussProd;
}
