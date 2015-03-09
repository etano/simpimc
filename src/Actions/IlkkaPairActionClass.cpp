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
  int nx_u, ny_u;
  in.Read("u/offDiag/nx", nx_u);
  in.Read("u/offDiag/ny", ny_u);
  vec<double> x_u(nx_u);
  vec<double> y_u(ny_u);
  mat<double> u_xy(nx_u,ny_u);
  in.Read("u/offDiag/x", x_u);
  in.Read("u/offDiag/y", y_u);
  in.Read("u/offDiag/u_xy", u_xy);

  // Spline u
  NUgrid* x_u_grid = create_general_grid(x_u.memptr(), x_u.size());
  NUgrid* y_u_grid = create_general_grid(y_u.memptr(), y_u.size());
  u_xy_spline = create_NUBspline_2d_d(x_u_grid, y_u_grid, xBC, xBC, u_xy.memptr());

  // u long range
  if (useLongRange) {
    // Read in r
    int nr_u;
    in.Read("u/diag/nrLong", nr_u);
    vec<double> r_u(nr_u);
    vec<double> uLong_r(nr_u);
    in.Read("u/diag/rLong", r_u);
    in.Read("u/diag/uLong_r", uLong_r);
    in.Read("u/diag/uLong_r0",uLong_r0);

    // Spline r
    NUgrid* r_u_grid = create_general_grid(r_u.memptr(), r_u.size());
    r_u_min = r_u_grid->start;
    r_u_max = r_u_grid->end;
    uLong_r_spline = create_NUBspline_1d_d(r_u_grid, xBC, uLong_r.memptr());

    // Read in k
    int nk_u;
    in.Read("u/diag/nk", nk_u);
    vec<double> k_u(nk_u);
    vec<double> tmpULong_k(nk_u);
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

  // Read in du
  int nx_du, ny_du;
  in.Read("du/offDiag/nx", nx_du);
  in.Read("du/offDiag/ny", ny_du);
  vec<double> x_du(nx_du);
  vec<double> y_du(ny_du);
  mat<double> du_xy(nx_du,ny_du);
  in.Read("du/offDiag/x", x_du);
  in.Read("du/offDiag/y", y_du);
  in.Read("du/offDiag/du_xy", du_xy);

  // Spline du
  NUgrid* x_du_grid = create_general_grid(x_du.memptr(), x_du.size());
  NUgrid* y_du_grid = create_general_grid(y_du.memptr(), y_du.size());
  du_xy_spline = create_NUBspline_2d_d(x_du_grid, y_du_grid, xBC, xBC, du_xy.memptr());

  // du long range
  if (useLongRange) {
    // Read in r
    int nr_du;
    in.Read("du/diag/nrLong", nr_du);
    vec<double> r_du(nr_du);
    vec<double> duLong_r(nr_du);
    in.Read("du/diag/rLong", r_du);
    in.Read("du/diag/duLong_r", duLong_r);
    in.Read("du/diag/duLong_r0",duLong_r0);

    // Spline r
    NUgrid* r_du_grid = create_general_grid(r_du.memptr(), r_du.size());
    r_du_min = r_du_grid->start;
    r_du_max = r_du_grid->end;
    duLong_r_spline = create_NUBspline_1d_d(r_du_grid, xBC, duLong_r.memptr());

    // Read in k
    int nk_du;
    in.Read("du/diag/nk", nk_du);
    vec<double> k_du(nk_du);
    vec<double> tmpDULong_k(nk_du);
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

  // Read in v
  int nr_v;
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
    int nr_vLong;
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
    int nk_v;
    in.Read("v/diag/nk", nk_v);
    vec<double> k_v(nk_v);
    vec<double> tmpVLong_k(nk_v);
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
    duLong_k0 *= 0.5*N1*N1*path.nBead;
    duLong_r0 *= -0.5*N1*path.nBead;
    vLong_k0 *= 0.5*N1*N1*path.nBead;
    vLong_r0 *= -0.5*N1*path.nBead;
  } else { // heterologous
    duLong_k0 *= N1*N2*path.nBead;
    duLong_r0 *= 0.*path.nBead;
    vLong_k0 *= N1*N2*path.nBead;
    vLong_r0 *= 0.*path.nBead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
double IlkkaPairAction::CalcV(double &r, double &rP, const int level)
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
double IlkkaPairAction::CalcU(double &r, double &rP, double &s, const int level)
{
  // Constants
  double q = 0.5*(r + rP);
  double x = q + 0.5*s;
  double y = q - 0.5*s;

  // Calculate U
  double U = 0.;
  eval_NUBspline_2d_d(u_xy_spline,x,y,&U);

  // Subtract out long range part
  if (useLongRange) {
    // Limits
    SetLimits(r_u_min, r_u_max, r, rP);

    // Splines
    double tmpU;
    eval_NUBspline_1d_d(uLong_r_spline,r,&tmpU);
    U -= 0.5*tmpU;
    eval_NUBspline_1d_d(uLong_r_spline,rP,&tmpU);
    U -= 0.5*tmpU;
  }

  return U;
}

/// Calculate the dU(r,r') value when given r and r' and the level 
double IlkkaPairAction::CalcdUdBeta(double &r, double &rP, double &s, const int level)
{
  // Constants
  double q = 0.5*(r + rP);
  double x = q + 0.5*s;
  double y = q - 0.5*s;

  // Calculate dU
  double dU = 0.;
  eval_NUBspline_2d_d(du_xy_spline,x,y,&dU);

  // Subtract out long range part
  if (useLongRange) {
    // Limits
    SetLimits(r_du_min, r_du_max, r, rP);

    // Splines
    double tmpDU;
    eval_NUBspline_1d_d(duLong_r_spline,r,&tmpDU);
    dU -= 0.5*tmpDU;
    eval_NUBspline_1d_d(duLong_r_spline,rP,&tmpDU);
    dU -= 0.5*tmpDU;
  }

  return dU;
}

/// Calculate the ULong value
double IlkkaPairAction::CalcVLong()
{
  // Get rho k
  field< vec< complex<double> > >& rhoK(path.GetRhoK());

  // Sum over k vectors
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iK=0; iK<path.ks.size(); iK++) {
    for (int iB=0; iB<path.nBead; iB++) {
      double rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
      tot += rhok2*vLong_k(iK);
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot + vLong_k0 + vLong_r0;
}

/// Calculate the ULong value
double IlkkaPairAction::CalcULong(const int b0, const int b1, const int level)
{
  // Get rho k
  field< vec< complex<double> > >& rhoK(path.GetRhoK());

  // Sum over k vectors
  int skip = 1<<level;
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iK=0; iK<path.ks.size(); iK++) {
    for (int iB=b0; iB<b1; iB+=skip) {
      double rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
      tot += uLong_k(iK)*rhok2;
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
double IlkkaPairAction::CalcdUdBetaLong()
{
  // Get rho k
  field< vec< complex<double> > >& rhoK(path.GetRhoK());

  // Sum over k vectors
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iK=0; iK<path.ks.size(); iK++) {
    for (int iB=0; iB<path.nBead; iB++) {
      double rhok2 = cmag2(rhoK(path.beadLoop(iB),iSpeciesA)(iK),rhoK(path.beadLoop(iB),iSpeciesB)(iK));
      tot += duLong_k(iK)*rhok2;
    }
  }

  if (iSpeciesB != iSpeciesA)
    tot *= 2.;

  return tot + duLong_k0 + duLong_r0;
}
