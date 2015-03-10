#include "DavidPairActionClass.h"

void DavidPairAction::ReadFile(string fileName)
{
  // Load file
  IOClass in;
  in.load(fileName);

  // Form strings
  stringstream tmpSS;
  tmpSS << "/Ukj" << nOrder;
  string UkjStr = tmpSS.str();
  stringstream tmpSS2;
  tmpSS2 << "/dUkjdBeta" << nOrder;
  string dUkjdBetaStr = tmpSS2.str();

  // Read in and create grid
  double rStart, rEnd;
  uint nGrid;
  string gridType;
  vec<double> gridPoints;
  in.Read(UkjStr + "/Grid/Start", rStart);
  in.Read(UkjStr + "/Grid/End", rEnd);
  in.Read(UkjStr + "/Grid/NumGridPoints", nGrid);
  in.Read(UkjStr + "/Grid/Type", gridType);
  if (gridType == "LOG")
    grid = create_log_grid(rStart, rEnd, nGrid);
  else if (gridType == "LINEAR")
    grid = create_linear_grid(rStart, rEnd, nGrid);
  else {
    in.Read(UkjStr + "/Grid/GridPoints", gridPoints);
    grid = create_general_grid(gridPoints.memptr(), nGrid);
  }

  // Read in taus
  nTau = maxLevel + 1;
  taus.set_size(nTau);
  in.Read(UkjStr + "/Taus", taus);
  bool tauFound = 0;
  for (uint iTau=0; iTau<nTau; iTau++)
    if (abs(taus(iTau)-path.tau) < 1.0e-6)
      tauFound = 1;
  if (!tauFound) {
    cerr << "ERROR: tau of " << path.tau << " not found." << endl;
    cerr << "Possible taus: " << taus << endl;
    exit(1);
  }

  // Read in potential
  vec<double> V(nGrid);
  in.Read("/Potential/Data", V);

  // Determine number of values for k,j sum
  nVal = 1;
  for (uint i = 1; i <= nOrder; ++i)
    nVal += 1+i;

  // Read in Ukj
  cube<double> tmpUkj(nVal,nGrid,nTau);
  in.Read(UkjStr + "/Data", tmpUkj);

  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Spline Ukj
  cube<double> tmpUkj2(nVal+1,nGrid,nTau);
  for(uint iTau=0; iTau<nTau; iTau++) {
    for (uint iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpUkj2(0,iGrid,iTau) = V(iGrid);
      for (uint iVal=1; iVal<nVal+1; ++iVal)
        tmpUkj2(iVal,iGrid,iTau) = tmpUkj(iVal-1,iGrid,iTau);
    }
    //for (uint iVal=1; iVal<nVal+1; ++iVal)
    //  tmpUkj2(iVal,nGrid-1,iTau) = 0.;
  }

  Ukj.set_size(nTau);
  for(uint iTau=0; iTau<nTau; iTau++) {
    Ukj(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    for (uint iVal=0; iVal<nVal+1; ++iVal) {
      vec<double> tmpV(nGrid);
      for (uint iGrid=0; iGrid<nGrid; ++iGrid)
        tmpV(iGrid) = tmpUkj2(iVal,iGrid,iTau);
      set_multi_NUBspline_1d_d(Ukj(iTau), iVal, tmpV.memptr());
    }
  }

  // Read in dUkjdBeta
  cube<double> tmpdUkjdBeta(nVal,nGrid,nTau);
  in.Read(dUkjdBetaStr + "/Data", tmpdUkjdBeta);

  // Spline dUkjdBeta
  cube<double> tmpdUkjdBeta2(nVal+1,nGrid,nTau);
  for(uint iTau=0; iTau<nTau; iTau++) {
    for (uint iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpdUkjdBeta2(0,iGrid,iTau) = V(iGrid);
      for (uint iVal=1; iVal<nVal+1; ++iVal)
        tmpdUkjdBeta2(iVal,iGrid,iTau) = tmpdUkjdBeta(iVal-1,iGrid,iTau);
    }
    //for (uint iVal=1; iVal<nVal+1; ++iVal)
    //  tmpdUkjdBeta2(iVal,nGrid-1,iTau) = 0.;
  }
  dUkjdBeta.set_size(nTau);
  for(uint iTau=0; iTau<nTau; iTau++) {
    dUkjdBeta(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    for (uint iVal=0; iVal<nVal+1; ++iVal) {
      vec<double> tmpV(nGrid);
      for (uint iGrid=0; iGrid<nGrid; ++iGrid)
        tmpV(iGrid) = tmpdUkjdBeta2(iVal,iGrid,iTau);
      set_multi_NUBspline_1d_d(dUkjdBeta(iTau), iVal, tmpV.memptr());
    }
  }

}

/// Calculate the U(r,r') value when given r and r' and the level 
double DavidPairAction::CalcV(double r, double rP, const uint level)
{
  // Limits
  double rMin, rMax;
  GetLimits(rMin, rMax, r, rP, grid);

  // This is the endpoint action
  double V;
  vec<double> rVals(nVal+1), rPVals(nVal+1), qVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rP,rPVals.memptr());
  V = 0.5*(rVals(0) + rPVals(0));

  return V;
}

/// Calculate the U(r,r') value when given r and r' and the level 
double DavidPairAction::CalcU(double r, double rP, double s, const uint level)
{
  // Constants
  double q = 0.5*(r + rP);
  double z = r - rP;

  // Limits
  double rMin, rMax;
  GetLimits(rMin, rMax, r, rP, grid);

  // This is the endpoint action
  double U;
  vec<double> rVals(nVal+1), rPVals(nVal+1), qVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rP,rPVals.memptr());
  U = 0.5*(rVals(1) + rPVals(1));

  // Add in off-diagonal terms
  if (s>0.0 && q<rMax) {
    vec<double> UqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    double z2 = z*z;
    double s2 = s*s;
    double s2inverse = 1./s2;
    double Sto2k = s2;
    for (uint k=1; k<=nOrder; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (uint j=0; j<=k; j++) {
        // indexing into the 2darray
        double Ucof = UqVals(k*(k+1)/2 + (j+1));
        U += (Ucof)*Zto2j*currS;
        Zto2j *= z2;
        currS *= s2inverse;
      }
      Sto2k *= s2;
    }
  }

  return U;

}

/// Calculate the U(r,r'), dU(r,r'), and V(r,r') value when given r and r' and the level 
double DavidPairAction::CalcdUdBeta(double r, double rP, double s, const uint level)
{
  // Constants
  double q = 0.5*(r + rP);
  double z = r - rP;

  // Limits
  double rMin, rMax;
  GetLimits(rMin, rMax, r, rP, grid);

  // This is the endpoint action
  double U, V, dU;
  vec<double> rVals(nVal+1), rPVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rP,rPVals.memptr());
  V = 0.5*(rVals(0) + rPVals(0));
  U = 0.5*(rVals(1) + rPVals(1));
  eval_multi_NUBspline_1d_d(dUkjdBeta(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(dUkjdBeta(level),rP,rPVals.memptr());
  dU = 0.5*(rVals(1) + rPVals(1));

  // Compensate for potential, which is subtracted from diaganal action in dm file.
  dU += V;

  // Add in off-diagonal terms
  if (s > 0.0 && q<rMax) {
    vec<double> UqVals(nVal+1), dUqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    eval_multi_NUBspline_1d_d(dUkjdBeta(level),q,dUqVals.memptr());
    double z2 = z*z;
    double s2 = s*s;
    double s2inverse = 1./s2;
    double Sto2k = s2;
    for (uint k=1; k<=nOrder; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (uint j=0; j<=k; j++){
        // indexing into the 2darray
        double Ucof = UqVals(k*(k+1)/2+j+1);
        double dUcof = dUqVals(k*(k+1)/2+j+1);
        U += (Ucof)*Zto2j*currS;
        dU += (dUcof)*Zto2j*currS;
        Zto2j *= z2;
        currS *= s2inverse;
      }
      Sto2k *= s2;
    }
  }

  return dU;
}
