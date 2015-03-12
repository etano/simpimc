#include "ilkka_pair_action_class.h"

void IlkkaPairAction::ReadFile(std::string file_name)
{
  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Load file
  IO in;
  in.Load(file_name);

  // Read in u
  uint n_x_u, n_y_u;
  in.Read("u/offDiag/n_x", n_x_u);
  in.Read("u/offDiag/n_y", n_y_u);
  vec<double> x_u(n_x_u);
  vec<double> y_u(n_y_u);
  mat<double> u_xy(n_x_u,n_y_u);
  in.Read("u/offDiag/x", x_u);
  in.Read("u/offDiag/y", y_u);
  in.Read("u/offDiag/u_xy", u_xy);

  // Spline u
  NUgrid* x_u_grid = create_general_grid(x_u.memptr(), x_u.size());
  NUgrid* y_u_grid = create_general_grid(y_u.memptr(), y_u.size());
  u_xy_spline = create_NUBspline_2d_d(x_u_grid, y_u_grid, xBC, xBC, u_xy.memptr());

  // u long range
  if (use_long_range) {
    // Read in r
    uint nr_u;
    in.Read("u/diag/nr_long", nr_u);
    vec<double> r_u(nr_u);
    vec<double> u_long_r(nr_u);
    in.Read("u/diag/r_long", r_u);
    in.Read("u/diag/u_long_r", u_long_r);
    in.Read("u/diag/u_long_r_0",u_long_r_0);

    // Spline r
    NUgrid* r_u_grid = create_general_grid(r_u.memptr(), r_u.size());
    r_u_min = r_u_grid->start;
    r_u_max = r_u_grid->end;
    u_long_r_spline = create_NUBspline_1d_d(r_u_grid, xBC, u_long_r.memptr());

    // Read in k
    uint n_k_u;
    in.Read("u/diag/n_k", n_k_u);
    vec<double> k_u(n_k_u);
    vec<double> tmp_u_long_k(n_k_u);
    in.Read("u/diag/k",k_u);
    in.Read("u/diag/u_long_k",tmp_u_long_k);
    in.Read("u/diag/u_long_k_0",u_long_k_0);

    // Build k std::vectors
    u_long_k.zeros(path.mag_ks.size());
    for (uint k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint k_i_u=0; k_i_u<k_u.size(); ++k_i_u) {
        if (fequal(path.mag_ks[k_i],k_u(k_i_u),1.e-8))
          u_long_k(k_i) = tmp_u_long_k(k_i_u);
      }
    }
  }

  // Read in du
  uint n_x_du, n_y_du;
  in.Read("du/offDiag/n_x", n_x_du);
  in.Read("du/offDiag/n_y", n_y_du);
  vec<double> x_du(n_x_du);
  vec<double> y_du(n_y_du);
  mat<double> du_xy(n_x_du,n_y_du);
  in.Read("du/offDiag/x", x_du);
  in.Read("du/offDiag/y", y_du);
  in.Read("du/offDiag/du_xy", du_xy);

  // Spline du
  NUgrid* x_du_grid = create_general_grid(x_du.memptr(), x_du.size());
  NUgrid* y_du_grid = create_general_grid(y_du.memptr(), y_du.size());
  du_xy_spline = create_NUBspline_2d_d(x_du_grid, y_du_grid, xBC, xBC, du_xy.memptr());

  // du long range
  if (use_long_range) {
    // Read in r
    uint nr_du;
    in.Read("du/diag/nr_long", nr_du);
    vec<double> r_du(nr_du);
    vec<double> du_long_r(nr_du);
    in.Read("du/diag/r_long", r_du);
    in.Read("du/diag/du_long_r", du_long_r);
    in.Read("du/diag/du_long_r_0",du_long_r_0);

    // Spline r
    NUgrid* r_du_grid = create_general_grid(r_du.memptr(), r_du.size());
    r_du_min = r_du_grid->start;
    r_du_max = r_du_grid->end;
    du_long_r_spline = create_NUBspline_1d_d(r_du_grid, xBC, du_long_r.memptr());

    // Read in k
    uint n_k_du;
    in.Read("du/diag/n_k", n_k_du);
    vec<double> k_du(n_k_du);
    vec<double> tmp_du_long_k(n_k_du);
    in.Read("du/diag/k",k_du);
    in.Read("du/diag/du_long_k",tmp_du_long_k);
    in.Read("du/diag/du_long_k_0",du_long_k_0);

    // Build k std::vectors
    du_long_k.zeros(path.mag_ks.size());
    for (uint k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint k_i_du=0; k_i_du<k_du.size(); ++k_i_du) {
        if (fequal(path.mag_ks[k_i],k_du(k_i_du),1.e-8))
          du_long_k(k_i) = tmp_du_long_k(k_i_du);
      }
    }

  }

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
  if (use_long_range) {
    // Read in r
    uint nr_v_long;
    in.Read("v/diag/nr_long", nr_v_long);
    vec<double> r_v_long(nr_v_long);
    vec<double> v_long_r(nr_v_long);
    in.Read("v/diag/r_long",r_v_long);
    in.Read("v/diag/v_long_r",v_long_r);
    in.Read("v/diag/v_long_r_0",v_long_r_0);

    // Spline r
    NUgrid* r_v_long_grid = create_general_grid(r_v_long.memptr(), r_v_long.size());
    r_v_long_min = r_v_long_grid->start;
    r_v_long_max = r_v_long_grid->end;
    v_long_r_spline = create_NUBspline_1d_d(r_v_long_grid, xBC, v_long_r.memptr());

    // Read in k
    uint n_k_v;
    in.Read("v/diag/n_k", n_k_v);
    vec<double> k_v(n_k_v);
    vec<double> tmp_v_long_k(n_k_v);
    in.Read("v/diag/k",k_v);
    in.Read("v/diag/v_long_k",tmp_v_long_k);
    in.Read("v/diag/v_long_k_0",v_long_k_0);

    // Build k std::vectors
    v_long_k.zeros(path.mag_ks.size());
    for (uint k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint k_i_v=0; k_i_v<k_v.size(); ++k_i_v) {
        if (fequal(path.mag_ks[k_i],k_v(k_i_v),1.e-8))
          v_long_k(k_i) = tmp_v_long_k(k_i_v);
      }
    }
  }

  // Calculate constants
  uint N1 = path.species_list[species_a_i]->n_part;
  uint N2 = path.species_list[species_b_i]->n_part;
  if (species_a_i == species_b_i) { // homologous
    du_long_k_0 *= 0.5*N1*N1*path.n_bead;
    du_long_r_0 *= -0.5*N1*path.n_bead;
    v_long_k_0 *= 0.5*N1*N1*path.n_bead;
    v_long_r_0 *= -0.5*N1*path.n_bead;
  } else { // heterologous
    du_long_k_0 *= N1*N2*path.n_bead;
    du_long_r_0 *= 0.*path.n_bead;
    v_long_k_0 *= N1*N2*path.n_bead;
    v_long_r_0 *= 0.*path.n_bead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
double IlkkaPairAction::CalcV(double r, double r_p, const uint level)
{
  // Limits
  SetLimits(r_v_min, r_v_max, r, r_p);

  // Calculate V
  double v = 0.;
  double tmp_v;
  eval_NUBspline_1d_d(v_r_spline,r,&tmp_v);
  v += 0.5*tmp_v;
  eval_NUBspline_1d_d(v_r_spline,r_p,&tmp_v);
  v += 0.5*tmp_v;
  if (use_long_range) {
    SetLimits(r_v_long_min, r_v_long_max, r, r_p);
    eval_NUBspline_1d_d(v_long_r_spline,r,&tmp_v);
    v -= 0.5*tmp_v;
    eval_NUBspline_1d_d(v_long_r_spline,r_p,&tmp_v);
    v -= 0.5*tmp_v;
  }

  return v;
}

/// Calculate the U(r,r') value when given r and r' and the level 
double IlkkaPairAction::CalcU(double r, double r_p, double s, const uint level)
{
  // Constants
  double q = 0.5*(r + r_p);
  double x = q + 0.5*s;
  double y = q - 0.5*s;

  // Calculate U
  double u = 0.;
  eval_NUBspline_2d_d(u_xy_spline,x,y,&u);

  // Subtract out long range part
  if (use_long_range) {
    // Limits
    SetLimits(r_u_min, r_u_max, r, r_p);

    // Splines
    double tmp_u;
    eval_NUBspline_1d_d(u_long_r_spline,r,&tmp_u);
    u -= 0.5*tmp_u;
    eval_NUBspline_1d_d(u_long_r_spline,r_p,&tmp_u);
    u -= 0.5*tmp_u;
  }

  return u;
}

/// Calculate the dU(r,r') value when given r and r' and the level 
double IlkkaPairAction::CalcdUdBeta(double r, double r_p, double s, const uint level)
{
  // Constants
  double q = 0.5*(r + r_p);
  double x = q + 0.5*s;
  double y = q - 0.5*s;

  // Calculate dU
  double du = 0.;
  eval_NUBspline_2d_d(du_xy_spline,x,y,&du);

  // Subtract out long range part
  if (use_long_range) {
    // Limits
    SetLimits(r_du_min, r_du_max, r, r_p);

    // Splines
    double tmp_du;
    eval_NUBspline_1d_d(du_long_r_spline,r,&tmp_du);
    du -= 0.5*tmp_du;
    eval_NUBspline_1d_d(du_long_r_spline,r_p,&tmp_du);
    du -= 0.5*tmp_du;
  }

  return du;
}

/// Calculate the ULong value
double IlkkaPairAction::CalcVLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint b_i=0; b_i<path.n_bead; b_i++) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += rhok2*v_long_k(k_i);
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + v_long_k_0 + v_long_r_0;
}

/// Calculate the ULong value
double IlkkaPairAction::CalcULong(const uint b0, const uint b1, const uint level)
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  uint skip = 1<<level;
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint b_i=b0; b_i<b1; b_i+=skip) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += u_long_k(k_i)*rhok2;
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
double IlkkaPairAction::CalcdUdBetaLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint b_i=0; b_i<path.n_bead; b_i++) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += du_long_k(k_i)*rhok2;
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + du_long_k_0 + du_long_r_0;
}