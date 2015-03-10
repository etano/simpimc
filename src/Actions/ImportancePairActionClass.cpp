#include "ImportancePairActionClass.h"

void ImportancePairAction::ReadFile(string fileName)
{
  // Set that this is an importance weight
  isImportanceWeight = 1;

  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Load file
  IOClass in;
  in.load(fileName);

    // Read in v
  uint nr_v;
  in.Read("v/diag/nr", nr_v);
  vec<double> r_v(nr_v);
  vec<double> v_r(nr_v);
  in.Read("v/diag/r", r_v);
  in.Read("v/diag/v_r", v_r);

  // Spline v
  NUgrid* r_v_grid = create_general_grid(r_v.memptr(), r_v.size());
  r_v_min = r_v_grid->start;
  r_v_max = r_v_grid->end;
  v_r_spline = create_NUBspline_1d_d(r_v_grid, xBC, v_r.memptr());

  // v long range
  if (useLongRange) {
    // Read in r
    uint nr_vLong;
    in.Read("v/diag/nrLong", nr_vLong);
    vec<double> r_vLong(nr_vLong);
    vec<double> vLong_r(nr_vLong);
    in.Read("v/diag/rLong",r_vLong);
    in.Read("v/diag/vLong_r",vLong_r);
    in.Read("v/diag/vLong_r0",vLong_r0);

    // Spline r
    NUgrid* r_vLong_grid = create_general_grid(r_vLong.memptr(), r_vLong.size());
    r_vLong_min = r_vLong_grid->start;
    r_vLong_max = r_vLong_grid->end;
    vLong_r_spline = create_NUBspline_1d_d(r_vLong_grid, xBC, vLong_r.memptr());

    // Read in k
    uint nk_v;
    in.Read("v/diag/nk", nk_v);
    vec<double> k_v(nk_v);
    vec<double> tmpVLong_k(nk_v);
    in.Read("v/diag/k",k_v);
    in.Read("v/diag/vLong_k",tmpVLong_k);
    in.Read("v/diag/vLong_k0",vLong_k0);

    // Build k vectors
    vLong_k.zeros(path.magKs.size());
    for (uint iK=0; iK<path.magKs.size(); ++iK) {
      for (uint iK_v=0; iK_v<k_v.size(); ++iK_v) {
        if (fequal(path.magKs[iK],k_v(iK_v),1.e-8))
          vLong_k(iK) = tmpVLong_k(iK_v);
      }
    }
  }

  // Calculate constants
  uint N1 = path.speciesList[iSpeciesA]->nPart;
  uint N2 = path.speciesList[iSpeciesB]->nPart;
  if (iSpeciesA == iSpeciesB) { // homologous
    vLong_k0 *= 0.5*N1*N1*path.nBead;
    vLong_r0 *= -0.5*N1*path.nBead;
  } else { // heterologous
    vLong_k0 *= N1*N2*path.nBead;
    vLong_r0 *= 0.*path.nBead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
double ImportancePairAction::CalcV(double r, double rP, const uint level)
{
  // Limits
  SetLimits(r_v_min, r_v_max, r, rP);

  // Calculate V
  double V = 0.;
  double tmpV;
  eval_NUBspline_1d_d(v_r_spline,r,&tmpV);
  V += 0.5*tmpV;
  eval_NUBspline_1d_d(v_r_spline,rP,&tmpV);
  V += 0.5*tmpV;
  if (useLongRange) {
    SetLimits(r_vLong_min, r_vLong_max, r, rP);
    eval_NUBspline_1d_d(vLong_r_spline,r,&tmpV);
    V -= 0.5*tmpV;
    eval_NUBspline_1d_d(vLong_r_spline,rP,&tmpV);
    V -= 0.5*tmpV;
  }

  return V;
}

/// Calculate the U(r,r') value when given r and r' and the level 
double ImportancePairAction::CalcU(double r, double rP, double s, const uint level)
{
  return CalcV(r,rP,level);
}

/// Calculate the dU(r,r') value when given r and r' and the level 
double ImportancePairAction::CalcdUdBeta(double r, double rP, double s, const uint level)
{
  return CalcV(r,rP,level);
}

/// Calculate the ULong value
double ImportancePairAction::CalcVLong()
{
  // Get rho k
  field< vec< complex<double> > >& rhoK(path.GetRhoK());

  // Sum over k vectors
  double tot = 0.;
  for (uint iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (uint iB=0; iB<path.nBead; iB++) {
        double rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += rhok2*vLong_k(iK);
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot + vLong_k0 + vLong_r0;
}

/// Calculate the ULong value
double ImportancePairAction::CalcULong(const uint b0, const uint b1, const uint level)
{
  // Get rho k
  field< vec< complex<double> > >& rhoK(path.GetRhoK());

  // Sum over k vectors
  uint skip = 1<<level;
  double tot = 0.;
  for (uint iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (uint iB=b0; iB<b1; iB+=skip) {
        double rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += vLong_k(iK)*rhok2;
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
double ImportancePairAction::CalcdUdBetaLong()
{
  return CalcVLong();
}
