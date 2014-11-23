#include "BarePairActionClass.h"

void BarePairAction::ReadFile(string fileName)
{
  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Load file
  IOClass in;
  in.load(fileName);

    // Read in v
  int nr_v;
  in.Read("v/diag/nr", nr_v);
  Tvector r_v(nr_v);
  Tvector v_r(nr_v);
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
    int nr_vLong;
    in.Read("v/diag/nrLong", nr_vLong);
    Tvector r_vLong(nr_vLong);
    Tvector vLong_r(nr_vLong);
    in.Read("v/diag/rLong",r_vLong);
    in.Read("v/diag/vLong_r",vLong_r);
    in.Read("v/diag/vLong_r0",vLong_r0);

    // Spline r
    NUgrid* r_vLong_grid = create_general_grid(r_vLong.memptr(), r_vLong.size());
    r_vLong_min = r_vLong_grid->start;
    r_vLong_max = r_vLong_grid->end;
    vLong_r_spline = create_NUBspline_1d_d(r_vLong_grid, xBC, vLong_r.memptr());

    // Read in k
    int nk_v;
    in.Read("v/diag/nk", nk_v);
    Tvector k_v(nk_v);
    Tvector tmpVLong_k(nk_v);
    in.Read("v/diag/k",k_v);
    in.Read("v/diag/vLong_k",tmpVLong_k);
    in.Read("v/diag/vLong_k0",vLong_k0);

    // Build k vectors
    vLong_k.zeros(path.magKs.size());
    for (int iK=0; iK<path.magKs.size(); ++iK) {
      for (int iK_v=0; iK_v<k_v.size(); ++iK_v) {
        if (fequal(path.magKs[iK],k_v(iK_v),1.e-8))
          vLong_k(iK) = tmpVLong_k(iK_v);
      }
    }
  }

  // Calculate constants
  int N1 = path.speciesList[iSpeciesA]->nPart;
  int N2 = path.speciesList[iSpeciesB]->nPart;
  if (iSpeciesA == iSpeciesB) { // homologous
    vLong_k0 *= 0.5*N1*N1*path.nBead;
    vLong_r0 *= -0.5*N1*path.nBead;
  } else { // heterologous
    vLong_k0 *= N1*N2*path.nBead;
    vLong_r0 *= 0.*path.nBead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
RealType BarePairAction::CalcV(RealType &r, RealType &rP, int level)
{
  // Limits
  SetLimits(r_v_min, r_v_max, r, rP);

  // Calculate V
  RealType V = 0.;
  RealType tmpV;
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
RealType BarePairAction::CalcU(RealType &r, RealType &rP, RealType &s, int level)
{
  int skip = 1>>level;
  RealType levelTau = skip*path.tau;
  return levelTau*CalcV(r,rP,level);
}

/// Calculate the dU(r,r') value when given r and r' and the level 
RealType BarePairAction::CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level)
{
  return CalcV(r,rP,level);
}

/// Calculate the ULong value
RealType BarePairAction::CalcVLong()
{
  // Get rho k
  field<Cvector>& rhoK(path.GetRhoK());

  // Sum over k vectors
  RealType tot = 0.;
  for (int iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (int iB=0; iB<path.nBead; iB++) {
        RealType rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += rhok2*vLong_k(iK);
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot + vLong_k0 + vLong_r0;
}

/// Calculate the ULong value
RealType BarePairAction::CalcULong(int b0, int b1, int level)
{
  // Get rho k
  field<Cvector>& rhoK(path.GetRhoK());

  // Sum over k vectors
  int skip = 1<<level;
  RealType tot = 0.;
  for (int iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (int iB=b0; iB<b1; iB+=skip) {
        RealType rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += vLong_k(iK)*rhok2;
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
RealType BarePairAction::CalcdUdBetaLong()
{
  return CalcVLong();
}
