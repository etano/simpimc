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
  RealType rStart, rEnd;
  int nGrid;
  string gridType;
  Tvector gridPoints;
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
  for (int iTau=0; iTau<nTau; iTau++)
    if (abs(taus(iTau)-path.tau) < 1.0e-6)
      tauFound = 1;
  if (!tauFound) {
    cerr << "ERROR: tau of " << path.tau << " not found." << endl;
    cerr << "Possible taus: " << taus << endl;
    exit(1);
  }

  // Read in potential
  Tvector V(nGrid);
  in.Read("/Potential/Data", V);

  // Determine number of values for k,j sum
  nVal = 1;
  for (int i = 1; i <= nOrder; ++i)
    nVal += 1+i;

  // Read in Ukj
  Tcube tmpUkj(nVal,nGrid,nTau);
  in.Read(UkjStr + "/Data", tmpUkj);

  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Spline Ukj
  Tcube tmpUkj2(nVal+1,nGrid,nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    for (int iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpUkj2(0,iGrid,iTau) = V(iGrid);
      for (int iVal=1; iVal<nVal+1; ++iVal)
        tmpUkj2(iVal,iGrid,iTau) = tmpUkj(iVal-1,iGrid,iTau);
    }
    //for (int iVal=1; iVal<nVal+1; ++iVal)
    //  tmpUkj2(iVal,nGrid-1,iTau) = 0.;
  }

  Ukj.set_size(nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    Ukj(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    for (int iVal=0; iVal<nVal+1; ++iVal) {
      Tvector tmpV(nGrid);
      for (int iGrid=0; iGrid<nGrid; ++iGrid)
        tmpV(iGrid) = tmpUkj2(iVal,iGrid,iTau);
      set_multi_NUBspline_1d_d(Ukj(iTau), iVal, tmpV.memptr());
    }
  }

  // Read in dUkjdBeta
  Tcube tmpdUkjdBeta(nVal,nGrid,nTau);
  in.Read(dUkjdBetaStr + "/Data", tmpdUkjdBeta);

  // Spline dUkjdBeta
  Tcube tmpdUkjdBeta2(nVal+1,nGrid,nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    for (int iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpdUkjdBeta2(0,iGrid,iTau) = V(iGrid);
      for (int iVal=1; iVal<nVal+1; ++iVal)
        tmpdUkjdBeta2(iVal,iGrid,iTau) = tmpdUkjdBeta(iVal-1,iGrid,iTau);
    }
    //for (int iVal=1; iVal<nVal+1; ++iVal)
    //  tmpdUkjdBeta2(iVal,nGrid-1,iTau) = 0.;
  }
  dUkjdBeta.set_size(nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    dUkjdBeta(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    for (int iVal=0; iVal<nVal+1; ++iVal) {
      Tvector tmpV(nGrid);
      for (int iGrid=0; iGrid<nGrid; ++iGrid)
        tmpV(iGrid) = tmpdUkjdBeta2(iVal,iGrid,iTau);
      set_multi_NUBspline_1d_d(dUkjdBeta(iTau), iVal, tmpV.memptr());
    }
  }

}

/// Calculate the U(r,r') value when given r and r' and the level 
RealType DavidPairAction::CalcU(RealType &r, RealType &rP, RealType &s, int level)
{
  // Constants
  RealType q = 0.5*(r + rP);
  RealType z = r - rP;

  // Limits
  RealType rMin, rMax;
  GetLimits(rMin, rMax, r, rP, grid);

  // This is the endpoint action
  RealType U;
  Tvector rVals(nVal+1), rPVals(nVal+1), qVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rP,rPVals.memptr());
  U = 0.5*(rVals(1) + rPVals(1));

  // Add in off-diagonal terms
  if (s>0.0 && q<rMax) {
    Tvector UqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    RealType z2 = z*z;
    RealType s2 = s*s;
    RealType s2inverse = 1./s2;
    RealType Sto2k = s2;
    for (int k=1; k<=nOrder; k++) {
      RealType Zto2j = 1;
      RealType currS = Sto2k;
      for (int j=0; j<=k; j++) {
        // indexing into the 2darray
        RealType Ucof = UqVals(k*(k+1)/2 + (j+1));
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
RealType DavidPairAction::CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level)
{
  // Constants
  RealType q = 0.5*(r + rP);
  RealType z = r - rP;

  // Limits
  RealType rMin, rMax;
  GetLimits(rMin, rMax, r, rP, grid);

  // This is the endpoint action
  RealType U, V, dU;
  Tvector rVals(nVal+1), rPVals(nVal+1);
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
    Tvector UqVals(nVal+1), dUqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    eval_multi_NUBspline_1d_d(dUkjdBeta(level),q,dUqVals.memptr());
    RealType z2 = z*z;
    RealType s2 = s*s;
    RealType s2inverse = 1./s2;
    RealType Sto2k = s2;
    for (int k=1; k<=nOrder; k++) {
      RealType Zto2j = 1;
      RealType currS = Sto2k;
      for (int j=0; j<=k; j++){
        // indexing into the 2darray
        RealType Ucof = UqVals(k*(k+1)/2+j+1);
        RealType dUcof = dUqVals(k*(k+1)/2+j+1);
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
