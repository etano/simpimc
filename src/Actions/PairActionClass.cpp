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

  // Read in grid
  RealType rStart, rEnd;
  int nGrid;
  string gridType;
  in.Read(UkjStr + "/Grid/Start", rStart);
  in.Read(UkjStr + "/Grid/Start", rEnd);
  in.Read(UkjStr + "/Grid/NumGridPoints", nGrid);
  in.Read(UkjStr + "/Grid/Type", gridType);

  // Create grid
  grid = create_log_grid(rStart, rEnd, nGrid);

  // Read in potential
  V.set_size(nGrid);
  in.Read("/Potential/Data", V);

  // Determine number of values for k,j sum
  int nVal = 1;
  for (int i=0; i<nOrder; ++i)
    nVal += 2+i;
  int nTau = maxLevel+1;

  // Read in Ukj
  Tcube tmpUkj(nGrid,nVal,nTau);
  in.Read(UkjStr + "/Data", tmpUkj);

  // Spline Ukj
  Tcube tmpUkj2(nGrid,nVal+1,nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    for (int iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpUkj2(iGrid,0,iTau) = V(iGrid);
      for (int iVal=1; iVal<nVal+1; ++iVal)
        tmpUkj2(iGrid,iVal,iTau) = tmpUkj(iGrid,iVal-1,iTau);
    }
    for (int iVal=1; iVal<nVal+1; ++iVal)
      tmpUkj2(nGrid-1,iVal,iTau) = 0.;
  }
  double startDerv(5.0e30), endDerv(0.);
  BCtype_d xBC = {DERIV1, FLAT, startDerv, endDerv};
  Ukj.set_size(nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    Ukj(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    set_multi_NUBspline_1d_d(Ukj(iTau), nVal+1, tmpUkj2.slice(iTau).memptr());
  }

  // Read in dUkjdBeta
  Tcube tmpdUkjdBeta(nGrid,nVal,maxLevel+1);
  in.Read(dUkjdBetaStr + "/Data", tmpdUkjdBeta);

  // Spline dUkjdBeta
  Tcube tmpdUkjdBeta2(nGrid,nVal+1,nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    for (int iGrid=0; iGrid<nGrid-1; ++iGrid) {
      tmpdUkjdBeta2(iGrid,0,iTau) = V(iGrid);
      for (int iVal=1; iVal<nVal+1; ++iVal)
        tmpdUkjdBeta2(iGrid,iVal,iTau) = tmpdUkjdBeta(iGrid,iVal-1,iTau);
    }
    for (int iVal=1; iVal<nVal+1; ++iVal)
      tmpdUkjdBeta2(nGrid-1,iVal,iTau) = 0.;
  }
  dUkjdBeta.set_size(nTau);
  for(int iTau=0; iTau<nTau; iTau++) {
    dUkjdBeta(iTau) = create_multi_NUBspline_1d_d(grid, xBC, nVal+1);
    set_multi_NUBspline_1d_d(dUkjdBeta(iTau), nVal+1, tmpdUkjdBeta2.slice(iTau).memptr());
  }

  // Read in taus
  taus.set_size(maxLevel+1);
  in.Read(UkjStr + "/Taus", taus);

}

RealType PairAction::DActionDBeta()
{
  RealType tot = 0.;
  Tvector dr;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
      for (int iB=0; iB<path.nBead; iB+=1) {
        path.Dr(path(iP,iB),path(jP,iB),dr);
        tot += 1/sqrt(dot(dr,dr));
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
        path.Dr(path(iP,iB),path(jP,iB),dr);
        tot += levelTau*1/sqrt(dot(dr,dr));
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

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::UdUVsqz(RealType s, RealType q, RealType z, int level, RealType &U, RealType &dU, RealType &V)
{
  RealType r = q + 0.5*z;
  RealType rPrime = q - 0.5*z;

  // Limits
  //RealType rMax = grid->End;
  //if (r > rMax)
  //  r = rMax;
  //if (rPrime > rMax)
  //  rPrime = rMax;
  //RealType rmin = grid->Start;
  //if (rPrime < rmin)
  //  rPrime = rmin;
  //if(r < rmin)
  //  r = rmin;

  // This is the endpoint action
  Tvector rVals(nOrder+2), rPrimeVals(nOrder+2);
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
  //if (s > 0.0 && q<rMax) {
  if (s > 0.0) {
    Tvector UqVals(nOrder+2), dUqVals(nOrder+2);
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
        RealType Ucof = tmpUkjArray(k*(k+1)/2+j+1);
        RealType dUcof = tmpdUkjdBetaArray(k*(k+1)/2+j+1);
        U += (Ucof)*Zto2j*currS;
        dU += (dUcof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void PairAction::Usqz(RealType s, RealType q, RealType z, int level, RealType &U)
{
  RealType r = q + 0.5*z;
  RealType rPrime = q - 0.5*z;

  // Limits
  //if (r > grid->End)
  //  r = grid->End;
  //RealType rMax = grid->End;
  //if (rPrime > rMax)
  //  rPrime = rMax;
  //RealType rMin = grid->Start;
  //if (rPrime < rMin)
  //  rPrime = rMin;
  //if(r < rMin)
  //  r = rMin;

  // This is the endpoint action
  Tvector rVals(nOrder+2), rPrimeVals(nOrder+2), qVals(nOrder+2);
  eval_multi_NUBspline_1d_d(Ukj(level),r,rVals.memptr());
  eval_multi_NUBspline_1d_d(Ukj(level),rPrime,rPrimeVals.memptr());
  U = 0.5*(rVals(1) + rPrimeVals(1));

  // Add in off-diagonal terms
  //if (s>0.0 && q<rMax) {
  if (s>0.0) {
    Tvector UqVals(nOrder+2);
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
        RealType Ucof = tmpUkjArray(k*(k+1)/2 + (j+1));
        U += (Ucof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}



