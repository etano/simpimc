#include "OptimizedNodalClass.h"

void OptimizedNodal::Init(Input &in)
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

  // Read in variational parameters
  iModel = in.getAttribute<int>("model");
  vector<Input> paramSetInputs = in.getChildList("ParamSet");
  for (auto& paramSetInput : paramSetInputs) {
    vector<Input> paramInputs = paramSetInput.getChildList("Param");
    vector<double> params;
    for (auto& paramInput : paramInputs)
      params.push_back(paramInput.getAttribute<double>("val"));
    paramSets.push_back(params);
  }

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

// Create spline
void OptimizedNodal::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  if (path.PBC) {
    r_grid.start = -path.L/2.;
    r_grid.end = path.L/2.;
  } else {
    r_grid.start = 100.;
    r_grid.end = 100.;
    nImages = 0;
  }
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  int nSpline = path.nBead/2 + (path.nBead%2) + 1;
  rho_node_r_splines.set_size(paramSets.size(), nSpline);

  // Create splines
  for (int iParamSet=0; iParamSet<paramSets.size(); ++iParamSet) {
    #pragma omp parallel for
    for (int iSpline=0; iSpline<nSpline; ++iSpline) {
      vec<double> rho_node_r(r_grid.num);
      double t_i4LambdaTau = Geti4LambdaTau(iSpline+1); // TODO: This is hard-coded for free-particle-like nodal structures.

      // Make rho_free
      for (int i=0; i<r_grid.num; ++i) {
        double r = r_grid.start + i*dr;
        double r2i4LambdaTau = r*r*t_i4LambdaTau;
        rho_node_r(i) = 0.;
        for (int image=-nImages; image<=nImages; image++) {
          if (image != 0) {
            double t_r = r + image*path.L;
            rho_node_r(i) += exp(r2i4LambdaTau - t_r*t_r*t_i4LambdaTau);
          }
        }
        rho_node_r(i) = log1p(min(10.,rho_node_r(i)));
      }
      BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
      UBspline_1d_d* rho_node_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_node_r.memptr());
      rho_node_r_splines(iParamSet,iSpline) = rho_node_r_spline;
    }
    cout << "...param set " << iParamSet << " complete." << endl;
  }

  // Reset iParamSet
  SetParamSet(0);
}

double OptimizedNodal::GetGij(const vec<double>& r, const int sliceDiff)
{
  double gaussProd = 1.;
  double t_i4LambdaTau = Geti4LambdaTau(sliceDiff);
  for (int iD=0; iD<path.nD; iD++) {
    double gaussSum;
    eval_UBspline_1d_d(rho_node_r_splines(iParamSet,sliceDiff-1),r(iD),&gaussSum);
    gaussSum = exp(0.9999*gaussSum);
    gaussSum *= exp(-(r(iD)*r(iD)*t_i4LambdaTau));
    gaussProd *= gaussSum;
  }
  return gaussProd;
}

double OptimizedNodal::Geti4LambdaTau(const int sliceDiff)
{
  double t_i4LambdaTau(i4LambdaTau);
  if (iModel == 0)
    t_i4LambdaTau *= paramSets[iParamSet][0]/sliceDiff;
  return t_i4LambdaTau;
}

void OptimizedNodal::Write() {}
