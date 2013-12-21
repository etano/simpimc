#include "EnergyClass.h"

void Energy::Init()
{
  Reset();
}

void Energy::Reset()
{
  nMeasure = 0;
  E = KE = VE = NE = 0.;
}

void Energy::Accumulate()
{
  KE += getKE(); // Add up total Kinetic energy
  VE += getVE(); // Add up total Potential energy
  NE += getNE(); // Add up total Nodal energy
  E = KE + VE + NE; // Add up total energy
  nMeasure += 1;
}

void Energy::Write()
{
  RealType norm = 1./(path.nBead*nMeasure);
  E *= norm;
  KE *= norm;
  VE *= norm;
  NE *= norm;

  if (FirstTime) {
    FirstTime = 0;
    out.CreateExtendableDataSet("/"+Name+"/", "Total", E);
    out.CreateExtendableDataSet("/"+Name+"/", "Kinetic", KE);
    out.CreateExtendableDataSet("/"+Name+"/", "Potential", VE);
    out.CreateExtendableDataSet("/"+Name+"/", "Nodal", NE);
  } else {
    out.AppendDataSet("/"+Name+"/", "Total", E);
    out.AppendDataSet("/"+Name+"/", "Kinetic", KE);
    out.AppendDataSet("/"+Name+"/", "Potential", VE);
    out.AppendDataSet("/"+Name+"/", "Nodal", NE);
  }

  Reset();
}

// Get Kinetic Energy
RealType Energy::getKE()
{
  int nImages = 0;

  RealType tot = 0.;
  Tvector dr(path.nD);
  for (int iP=0; iP<path.nPart; iP+=1) {
    RealType i4LambdaTau = 1./(4.*path.bead(iP,0)->species.lambda*path.tau);
    for (int iB=0; iB<path.nBead; iB+=1) {
      dr = path.bead(iP,iB)->r - path.bead(iP,iB)->next->r;
      path.PutInBox(dr);
      RealType gaussProd = 1.;
      Tvector gaussSum, numSum;
      gaussSum.zeros(path.nD);
      numSum.zeros(path.nD);
      for (int iD=0; iD<path.nD; iD++) {
        for (int image=-nImages; image<=nImages; image++) {
          RealType dist = dr(iD) + (RealType)image*path.L;
          RealType dist2i4LambdaTau = dist*dist*i4LambdaTau;
          RealType expPart = exp(-dist2i4LambdaTau);
          gaussSum(iD) += expPart;
          numSum(iD) += -(dist2i4LambdaTau/path.tau)*expPart;
        }
        gaussProd *= gaussSum(iD);
      }
      RealType scalarNumSum = 0.;
      for (int iD=0; iD<path.nD; iD++) {
        Tvector numProd;
        numProd.ones(path.nD);
        for (int jD=0; jD<path.nD; jD++) {
          if (iD != jD)
            numProd(iD) *= gaussSum(jD);
          else
            numProd(iD) *= numSum(jD);
        }
        scalarNumSum += numProd(iD);
      }
      tot += scalarNumSum/gaussProd;
    }
  }

  return path.nPartnBeadnDOver2Tau + tot;

}

// Get Potential Energy
inline RealType Energy::getVE()
{
  return getVEext() + getVEint();
}

// Get Interaction Potential Energy
RealType Energy::getVEint()
{
  RealType tot = 0.0;

  for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
    for (unsigned int iPart = 0; iPart < path.nPart-1; iPart += 1)  {
      for (unsigned int jPart = iPart+1; jPart < path.nPart; jPart += 1) {
        tot += 2.0 * path.getVint( path.bead(iPart,iBead) , path.bead(jPart,iBead) );
      }
    }
  }

  return tot;
}

// Get External Potential Energy
RealType Energy::getVEext()
{
  RealType tot = 0.0;

  if(path.trap) {
    for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
      for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
        tot += dot( path.bead(iPart,iBead) -> r , path.bead(iPart,iBead) -> r );
      }
    }
  }

  return path.halfOmega2 * path.onePlus3Tau2Omega2Over12 * tot;
}

// Get Nodal Energy
RealType Energy::getNE()
{
  if(!path.useNodeDist) return 0;

  //path.mode = 1;
  //for (unsigned int iPart = 0; iPart < path.nPart; iPart++) {
  //  for (unsigned int iBead = 0; iBead < path.nBead; iBead++) {
  //    path.updateNodeDistance(path.bead(iPart,iBead));
  //  }
  //}

  RealType xi, tot;
  tot = 0.0;
  for (unsigned int iPart = 0; iPart < path.nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < path.nBead; iBead += 1)  {
      RealType nD1 = path.bead(iPart,iBead) -> nDist;
      RealType nD2 = path.bead(iPart,iBead) -> next -> nDist;
      RealType nD1nD2 = nD1 * nD2;
      if (!nD1) {
        nD1nD2 = nD2 * nD2;
      } else if (!nD2) {
        nD1nD2 = nD1 * nD1;
      }
      //xi = nD1nD2*path.iLamTau;
      xi = 0.5*nD1nD2/(path.bead(iPart,iBead)->species.lambda*path.tau);
      tot += xi/(path.tau*expm1(xi));
      if (std::isnan(tot))
        std::cerr << tot << endl;
    }
  }

  return tot;

}
