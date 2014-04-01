#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder");
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  maxLevel = in.getAttribute<int>("maxLevel");
  GetOffset(speciesA,iSpeciesA,offsetA);
  GetOffset(speciesB,iSpeciesB,offsetB);

  string fileName = in.getAttribute<string>("file");
  ReadFile(fileName);

  out.Write("/Actions/"+name+"/file", fileName);
  out.Write("/Actions/"+name+"/nImages", nImages);
  out.Write("/Actions/"+name+"/nOrder", nOrder);
  out.Write("/Actions/"+name+"/speciesA", speciesA);
  out.Write("/Actions/"+name+"/speciesB", speciesB);
}

void PairAction::ReadFile(string fileName)
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
  in.Read(UkjStr + "/Grid/Start", rStart);
  in.Read(UkjStr + "/Grid/End", rEnd);
  in.Read(UkjStr + "/Grid/NumGridPoints", nGrid);
  in.Read(UkjStr + "/Grid/Type", gridType);
  grid = create_log_grid(rStart, rEnd, nGrid);

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

  RealType U, dU, v;
  RealType x(0.0), end(3.0);
  while (x<=end) {
    CalcUdUVsqz(0,x,0,0,U,dU,v);
    cout << x << " " << U << " " << dU << " " << v << endl;
    x += 0.1;
  }
  Tvector r;
  r.zeros(path.nD);
  Tvector rP;
  rP.zeros(path.nD);
  rP(0) = 0.5;
  CalcUdUVr(r,rP,0,U,dU,v);
  cout << " " << U << " " << dU << " " << v << endl;

}

RealType PairAction::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
      for (int iB=0; iB<path.nBead; iB+=1) {
        int jB = iB + 1;
        Tvector r, rP;
        path.Dr(path(iP,iB),path(jP,iB),r);
        path.Dr(path(iP,jB),path(jP,jB),rP);
        RealType U, dU, V;
        CalcUdUVr(r,rP,0,U,dU,V);
        tot += dU;
        //RealType sum = 0.;
        //for (int iD=0; iD<path.nD; iD++) {
        //  for (int image=-nImages; image<=nImages; image++) {
        //    RealType dist = dr(iD) + (RealType)image*path.L;
        //    sum += -0.5*1/dist;
        //  }
        //}
        //tot -= sum;
      }
    }
  }

  return tot;
}

RealType PairAction::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  if (level > maxLevel)
    return 0.;
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  Tvector dr;
  for (int iP=0; iP<particles.size(); ++iP) {
    for (int iB=b0; iB<b1; iB+=skip) {
      int jB = iB + skip;
      int iS = 0;
      int offset = 0;
      if (path(iP,iB)->species.name == speciesA) {
        iS = iSpeciesB;
        offset = offsetB;
      } else if (path(iP,iB)->species.name == speciesB) {
        iS = iSpeciesA;
        offset = offsetA;
      } else
        cerr << "ERROR: Unrecognized species in PairAction action." << endl;
      for (int jP=offset; jP<offset+path.speciesList[iS]->nPart; ++jP) {
        Tvector r, rP;
        path.Dr(path(iP,iB),path(jP,iB),r);
        path.Dr(path(iP,jB),path(jP,jB),rP);
        RealType U;
        CalcUr(r,rP,level,U);
        tot += U;
        //RealType gaussProd = 1.;
        //for (int iD=0; iD<path.nD; iD++) {
        //  RealType gaussSum = 0.;
        //  for (int image=-nImages; image<=nImages; image++) {
        //    RealType dist = dr(iD) + (RealType)image*path.L;
        //    gaussSum += exp(-levelTau*1/dist);
        //  }
        //  gaussProd *= gaussSum;
        //}
        //tot -= log(gaussProd);
      }
    }
  }

  return tot;
}

void PairAction::Write()
{

}

/// Calculate the U(r,r') value when given r and r' and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::CalcUdUVr(Tvector& r, Tvector& rP, int level, RealType &U, RealType &dU, RealType &V)
{
  Tvector tmp;
  RealType rMag = mag(r);
  RealType rPMag = mag(rP);

  RealType q = (rMag + rPMag)/2.;
  Tvector dr;
  path.Dr(r,rP,dr);
  RealType s = mag(dr);
  RealType z = rMag - rPMag;

  return CalcUdUVsqz(s, q, z, level, U, dU, V);
}

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::CalcUdUVsqz(RealType s, RealType q, RealType z, int level, RealType &U, RealType &dU, RealType &V)
{
  RealType r = q + 0.5*z;
  RealType rPrime = q - 0.5*z;

  // Limits
  RealType rMax = grid->end;
  if (r > rMax)
    r = rMax;
  if (rPrime > rMax)
    rPrime = rMax;
  RealType rmin = grid->start;
  if (rPrime < rmin)
    rPrime = rmin;
  if(r < rmin)
    r = rmin;

  // This is the endpoint action
  Tvector rVals(nVal+1), rPrimeVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rPrime,rPrimeVals.memptr());
  V = 0.5*(rVals(0) + rPrimeVals(0));
  U = 0.5*(rVals(1) + rPrimeVals(1));
  eval_multi_NUBspline_1d_d(dUkjdBeta(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(dUkjdBeta(level),rPrime,rPrimeVals.memptr());
  dU = 0.5*(rVals(1) + rPrimeVals(1));

  // Compensate for potential, which is subtracted from diaganal action in dm file.
  dU += V;

  // Add in off-diagonal terms
  if (s > 0.0 && q<rMax) {
    Tvector UqVals(nVal+1), dUqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    eval_multi_NUBspline_1d_d(dUkjdBeta(level),q,dUqVals.memptr());
    RealType zsquared = z*z;
    RealType ssquared = s*s;
    RealType ssquaredinverse = 1./ssquared;
    RealType Sto2k = ssquared;
    for (int k=1; k<=nOrder; k++) {
      RealType Zto2j = 1;
      RealType currS = Sto2k;
      for (int j=0; j<=k; j++){
        // indexing into the 2darray
        RealType Ucof = UqVals(k*(k+1)/2+j+1);
        RealType dUcof = dUqVals(k*(k+1)/2+j+1);
        U += (Ucof)*Zto2j*currS;
        dU += (dUcof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}

/// Calculate the U(r,r') value when given r and r' and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::CalcUr(Tvector& r, Tvector& rP, int level, RealType &U)
{
  Tvector tmp;
  RealType rMag = mag(r);
  RealType rPMag = mag(rP);

  RealType q = (rMag + rPMag)/2.;
  Tvector dr = r-rP;
  RealType s = mag(dr);
  RealType z = rMag - rPMag;
  CalcUsqz(s, q, z, level, U);
}

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::CalcUsqz(RealType s, RealType q, RealType z, int level, RealType &U)
{
  RealType r = q + 0.5*z;
  RealType rPrime = q - 0.5*z;

  // Limits
  RealType rMax = grid->end;
  if (r > rMax)
    r = rMax;
  if (rPrime > rMax)
    rPrime = rMax;
  RealType rMin = grid->start;
  if (rPrime < rMin)
    rPrime = rMin;
  if(r < rMin)
    r = rMin;

  // This is the endpoint action
  Tvector rVals(nVal+1), rPrimeVals(nVal+1), qVals(nVal+1);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rPrime,rPrimeVals.memptr());
  U = 0.5*(rVals(1) + rPrimeVals(1));

  // Add in off-diagonal terms
  if (s>0.0 && q<rMax) {
    Tvector UqVals(nVal+1);
    eval_multi_NUBspline_1d_d(Ukj(level),q,UqVals.memptr());
    RealType zsquared = z*z;
    RealType ssquared = s*s;
    RealType ssquaredinverse = 1./ssquared;
    RealType Sto2k = ssquared;
    for (int k=1; k<=nOrder; k++) {
      RealType Zto2j = 1;
      RealType currS = Sto2k;
      for (int j=0; j<=k; j++) {
        // indexing into the 2darray
        RealType Ucof = UqVals(k*(k+1)/2 + (j+1));
        U += (Ucof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}

