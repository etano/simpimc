#include "IlkkaPairActionClass.h"

void IlkkaPairAction::ReadFile(string fileName)
{
  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Load file
  IOClass in;
  in.load(fileName);

  nOrder = -1; // HACK

  // Read in u
  if (useLongRange) {
    int nr_u;
    in.Read("u/diag/nrLong", nr_u);
    r_u.set_size(nr_u);
    uLong_r.set_size(nr_u);
    in.Read("u/diag/rLong", r_u);
    in.Read("u/diag/uLong_r", uLong_r);
    in.Read("u/diag/uLong_r0",uLong_r0);
    int nk_u;
    in.Read("u/diag/nk", nk_u);
    Tvector k_u(nk_u);
    Tvector tmpULong_k(nk_u);
    in.Read("u/diag/k",k_u);
    in.Read("u/diag/uLong_k",tmpULong_k);
    in.Read("u/diag/uLong_k0",uLong_k0);

    // Build k vectors
    uLong_k.zeros(path.magKs.size());
    for (int iK=0; iK<path.magKs.size(); ++iK) {
      for (int iK_u=0; iK_u<k_u.size(); ++iK_u) {
        if (fequal(path.magKs[iK],k_u(iK_u),1.e-8))
          uLong_k(iK) = tmpULong_k(iK_u);
      }
    }
  }
  int nx_u, ny_u;
  in.Read("u/offDiag/nx", nx_u);
  in.Read("u/offDiag/ny", ny_u);
  x_u.set_size(nx_u);
  y_u.set_size(ny_u);
  u_xy.set_size(nx_u,ny_u);
  in.Read("u/offDiag/x", x_u);
  in.Read("u/offDiag/y", y_u);
  in.Read("u/offDiag/u_xy", u_xy);

  // Spline u
  x_u_grid = create_general_grid(x_u.memptr(), x_u.size());
  y_u_grid = create_general_grid(y_u.memptr(), y_u.size());
  u_xy_spline = create_NUBspline_2d_d(x_u_grid, y_u_grid, xBC, xBC, u_xy.memptr());

  if (useLongRange) {
    r_u_grid = create_general_grid(r_u.memptr(), r_u.size());
    uLong_r_spline = create_NUBspline_1d_d(r_u_grid, xBC, uLong_r.memptr());
  }

  // Read in du
  if (useLongRange) {
    int nr_du;
    in.Read("du/diag/nrLong", nr_du);
    r_du.set_size(nr_du);
    duLong_r.set_size(nr_du);
    in.Read("du/diag/rLong", r_du);
    in.Read("du/diag/duLong_r", duLong_r);
    in.Read("du/diag/duLong_r0",duLong_r0);
    int nk_du;
    in.Read("du/diag/nk", nk_du);
    Tvector k_du(nk_du);
    Tvector tmpDULong_k(nk_du);
    in.Read("du/diag/k",k_du);
    in.Read("du/diag/duLong_k",tmpDULong_k);
    in.Read("du/diag/duLong_k0",duLong_k0);

    // Build k vectors
    duLong_k.zeros(path.magKs.size());
    for (int iK=0; iK<path.magKs.size(); ++iK) {
      for (int iK_du=0; iK_du<k_du.size(); ++iK_du) {
        if (fequal(path.magKs[iK],k_du(iK_du),1.e-8))
          duLong_k(iK) = tmpDULong_k(iK_du);
      }
    }
  }
  int nx_du, ny_du;
  in.Read("du/offDiag/nx", nx_du);
  in.Read("du/offDiag/ny", ny_du);
  x_du.set_size(nx_du);
  y_du.set_size(ny_du);
  du_xy.set_size(nx_du,ny_du);
  in.Read("du/offDiag/x", x_du);
  in.Read("du/offDiag/y", y_du);
  in.Read("du/offDiag/du_xy", du_xy);

  // Spline du
  x_du_grid = create_general_grid(x_du.memptr(), x_du.size());
  y_du_grid = create_general_grid(y_du.memptr(), y_du.size());
  du_xy_spline = create_NUBspline_2d_d(x_du_grid, y_du_grid, xBC, xBC, du_xy.memptr());
  if (useLongRange) {
    r_du_grid = create_general_grid(r_du.memptr(), r_du.size());
    duLong_r_spline = create_NUBspline_1d_d(r_du_grid, xBC, duLong_r.memptr());
  }

  // Read in v
  int nr_v;
  in.Read("v/diag/nr", nr_v);
  r_v.set_size(nr_v);
  v_r.set_size(nr_v);
  in.Read("v/diag/r", r_v);
  in.Read("v/diag/v_r", v_r);
  if (useLongRange) {
    int nr_vLong;
    in.Read("v/diag/nrLong", nr_vLong);
    r_vLong.set_size(nr_vLong);
    vLong_r.set_size(nr_vLong);
    in.Read("v/diag/rLong",r_vLong);
    in.Read("v/diag/vLong_r",vLong_r);
    in.Read("v/diag/vLong_r0",vLong_r0);
    int nk_v;
    in.Read("du/diag/nk", nk_v);
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

  // Spline v
  r_v_grid = create_general_grid(r_v.memptr(), r_v.size());
  v_r_spline = create_NUBspline_1d_d(r_v_grid, xBC, v_r.memptr());
  if (useLongRange) {
    r_vLong_grid = create_general_grid(r_vLong.memptr(), r_vLong.size());
    vLong_r_spline = create_NUBspline_1d_d(r_vLong_grid, xBC, vLong_r.memptr());
  }

  // Calculate constants
  int N1 = path.speciesList[iSpeciesA]->nPart;
  int N2 = path.speciesList[iSpeciesB]->nPart;
  if (iSpeciesA == iSpeciesB) { // homologous
    duLong_k0 *= 0.25*N1*N1*path.nBead;
    duLong_r0 *= -0.5*N1*path.nBead;
    vLong_k0 *= 0.25*N1*N1*path.nBead;
    vLong_r0 *= -0.5*N1*path.nBead;
  } else { // heterologous
    duLong_k0 *= 0.5*N1*N2*path.nBead;
    duLong_r0 *= 0.*path.nBead;
    vLong_k0 *= 0.5*N1*N2*path.nBead;
    vLong_r0 *= 0.*path.nBead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
RealType IlkkaPairAction::CalcV(Tvector& rVec, Tvector& rPVec, int level)
{
  // Constants
  RealType r, rP, q, s, z, x, y;
  GetConstants(rVec, rPVec, r, rP, q, s, z, x, y);

  // Limits
  RealType rMin, rMax;
  GetLimits(rMin, rMax, r, rP, r_v_grid);

  // Calculate V
  RealType V = 0.;
  RealType tmpV;
  eval_NUBspline_1d_d(v_r_spline,r,&tmpV);
  V += 0.5*tmpV;
  eval_NUBspline_1d_d(v_r_spline,rP,&tmpV);
  V += 0.5*tmpV;
  if (useLongRange) {
    GetLimits(rMin, rMax, r, rP, r_vLong_grid);
    eval_NUBspline_1d_d(vLong_r_spline,r,&tmpV);
    V -= 0.5*tmpV;
    eval_NUBspline_1d_d(vLong_r_spline,rP,&tmpV);
    V -= 0.5*tmpV;
  }

  return V;
}

/// Calculate the U(r,r') value when given r and r' and the level 
RealType IlkkaPairAction::CalcU(RealType &r, RealType &rP, RealType &s, int level)
{
  // Constants
  RealType q = 0.5*(r + rP);
  RealType x = q + 0.5*s;
  RealType y = q - 0.5*s;

  // Calculate U
  RealType U = 0.;
  eval_NUBspline_2d_d(u_xy_spline,x,y,&U);

  // Subtract out long range part
  if (useLongRange) {
    // Limits
    RealType rMin, rMax;
    GetLimits(rMin, rMax, r, rP, r_u_grid);

    // Splines
    RealType tmpU;
    eval_NUBspline_1d_d(uLong_r_spline,r,&tmpU);
    U -= 0.5*tmpU;
    eval_NUBspline_1d_d(uLong_r_spline,rP,&tmpU);
    U -= 0.5*tmpU;
  }

  return U;
}

/// Calculate the dU(r,r') value when given r and r' and the level 
RealType IlkkaPairAction::CalcdUdBeta(RealType &r, RealType &rP, RealType &s, int level)
{
  // Constants
  RealType q = 0.5*(r + rP);
  RealType x = q + 0.5*s;
  RealType y = q - 0.5*s;

  // Calculate dU
  RealType dU = 0.;
  eval_NUBspline_2d_d(du_xy_spline,x,y,&dU);

  // Subtract out long range part
  if (useLongRange) {
    // Limits
    RealType rMin, rMax;
    GetLimits(rMin, rMax, r, rP, r_du_grid);

    // Splines
    RealType tmpDU;
    eval_NUBspline_1d_d(duLong_r_spline,r,&tmpDU);
    dU -= 0.5*tmpDU;
    eval_NUBspline_1d_d(duLong_r_spline,rP,&tmpDU);
    dU -= 0.5*tmpDU;
  }

  return dU;
}

/// Calculate the ULong value
RealType IlkkaPairAction::CalcVLong()
{
  // Get rho k
  arma::field<Cvector>& rhoK(path.GetRhoK());

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
RealType IlkkaPairAction::CalcULong(int b0, int b1, vector<int> &particles, int level)
{
  // Get rho k
  arma::field<Cvector>& rhoK(path.GetRhoK());

  // Sum over k vectors
  int skip = 1<<level;
  RealType tot = 0.;
  for (int iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (int iB=b0; iB<b1; iB+=skip) {
        RealType rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += rhok2*uLong_k(iK);
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
RealType IlkkaPairAction::CalcdUdBetaLong()
{
  // Get rho k
  arma::field<Cvector>& rhoK(path.GetRhoK());

  // Sum over k vectors
  RealType tot = 0.;
  for (int iK=0; iK<path.ks.size(); iK++) {
    if (path.magKs[iK] < kCut) {
      for (int iB=0; iB<path.nBead; iB++) {
        RealType rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
        tot += duLong_k(iK)*rhok2;
      }
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot + duLong_k0 + duLong_r0;
}
