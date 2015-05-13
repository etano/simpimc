#include "ilkka_pair_action_class.h"

void IlkkaPairAction::ReadFile(std::string file_name)
{
  // Boundary conditions
  BCtype_d xBC = {NATURAL, NATURAL}; // HACK: Is this correct?

  // Load file
  IO in;
  in.Load(file_name);

  // Read in u
  uint32_t n_x_u, n_y_u;
  in.Read("u/off_diag/n_x", n_x_u);
  in.Read("u/off_diag/n_y", n_y_u);
  vec<double> x_u(n_x_u);
  vec<double> y_u(n_y_u);
  mat<double> u_xy(n_x_u,n_y_u);
  in.Read("u/off_diag/x", x_u);
  in.Read("u/off_diag/y", y_u);
  in.Read("u/off_diag/u_xy", u_xy);

  // Spline u
  NUgrid* x_u_grid = create_general_grid(x_u.memptr(), x_u.size());
  NUgrid* y_u_grid = create_general_grid(y_u.memptr(), y_u.size());
  u_xy_spline = create_NUBspline_2d_d(x_u_grid, y_u_grid, xBC, xBC, u_xy.memptr());

  // u long range
  if (use_long_range) {
    // Read in r
    uint32_t n_r_u;
    in.Read("u/diag/n_r_long", n_r_u);
    vec<double> r_u(n_r_u);
    vec<double> u_long_r(n_r_u);
    in.Read("u/diag/r_long", r_u);
    in.Read("u/diag/u_long_r", u_long_r);
    in.Read("u/diag/u_long_r_0",u_long_r_0);

    // Spline r
    NUgrid* r_u_grid = create_general_grid(r_u.memptr(), r_u.size());
    r_u_min = r_u_grid->start;
    r_u_max = r_u_grid->end;
    u_long_r_spline = create_NUBspline_1d_d(r_u_grid, xBC, u_long_r.memptr());

    // Read in k
    uint32_t n_k_u;
    in.Read("u/diag/n_k", n_k_u);
    vec<double> k_u(n_k_u);
    vec<double> tmp_u_long_k(n_k_u);
    in.Read("u/diag/k",k_u);
    in.Read("u/diag/u_long_k",tmp_u_long_k);
    in.Read("u/diag/u_long_k_0",u_long_k_0);

    // Build k std::vectors
    u_long_k.zeros(path.mag_ks.size());
    for (uint32_t k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint32_t k_i_u=0; k_i_u<k_u.size(); ++k_i_u) {
        if (fequal(path.mag_ks[k_i],k_u(k_i_u),1.e-8))
          u_long_k(k_i) = tmp_u_long_k(k_i_u);
      }
    }
  }

  // Read in du
  uint32_t n_x_du, n_y_du;
  in.Read("du/off_diag/n_x", n_x_du);
  in.Read("du/off_diag/n_y", n_y_du);
  vec<double> x_du(n_x_du);
  vec<double> y_du(n_y_du);
  mat<double> du_xy(n_x_du,n_y_du);
  in.Read("du/off_diag/x", x_du);
  in.Read("du/off_diag/y", y_du);
  in.Read("du/off_diag/du_xy", du_xy);

  // Spline du
  NUgrid* x_du_grid = create_general_grid(x_du.memptr(), x_du.size());
  NUgrid* y_du_grid = create_general_grid(y_du.memptr(), y_du.size());
  du_xy_spline = create_NUBspline_2d_d(x_du_grid, y_du_grid, xBC, xBC, du_xy.memptr());

  // du long range
  if (use_long_range) {
    // Read in r
    uint32_t n_r_du;
    in.Read("du/diag/n_r_long", n_r_du);
    vec<double> r_du(n_r_du);
    vec<double> du_long_r(n_r_du);
    in.Read("du/diag/r_long", r_du);
    in.Read("du/diag/du_long_r", du_long_r);
    in.Read("du/diag/du_long_r_0",du_long_r_0);

    // Spline r
    NUgrid* r_du_grid = create_general_grid(r_du.memptr(), r_du.size());
    r_du_min = r_du_grid->start;
    r_du_max = r_du_grid->end;
    du_long_r_spline = create_NUBspline_1d_d(r_du_grid, xBC, du_long_r.memptr());

    // Read in k
    uint32_t n_k_du;
    in.Read("du/diag/n_k", n_k_du);
    vec<double> k_du(n_k_du);
    vec<double> tmp_du_long_k(n_k_du);
    in.Read("du/diag/k",k_du);
    in.Read("du/diag/du_long_k",tmp_du_long_k);
    in.Read("du/diag/du_long_k_0",du_long_k_0);

    // Build k std::vectors
    du_long_k.zeros(path.mag_ks.size());
    for (uint32_t k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint32_t k_i_du=0; k_i_du<k_du.size(); ++k_i_du) {
        if (fequal(path.mag_ks[k_i],k_du(k_i_du),1.e-8))
          du_long_k(k_i) = tmp_du_long_k(k_i_du);
      }
    }

  }

  // Read in v
  uint32_t n_r_v;
  in.Read("v/diag/n_r", n_r_v);
  vec<double> r_v(n_r_v);
  vec<double> v_r(n_r_v);
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
    uint32_t n_r_v_long;
    in.Read("v/diag/n_r_long", n_r_v_long);
    vec<double> r_v_long(n_r_v_long);
    vec<double> v_long_r(n_r_v_long);
    in.Read("v/diag/r_long",r_v_long);
    in.Read("v/diag/v_long_r",v_long_r);
    in.Read("v/diag/v_long_r_0",v_long_r_0);

    // Spline r
    NUgrid* r_v_long_grid = create_general_grid(r_v_long.memptr(), r_v_long.size());
    r_v_long_min = r_v_long_grid->start;
    r_v_long_max = r_v_long_grid->end;
    v_long_r_spline = create_NUBspline_1d_d(r_v_long_grid, xBC, v_long_r.memptr());

    // Read in k
    uint32_t n_k_v;
    in.Read("v/diag/n_k", n_k_v);
    vec<double> k_v(n_k_v);
    vec<double> tmp_v_long_k(n_k_v);
    in.Read("v/diag/k",k_v);
    in.Read("v/diag/v_long_k",tmp_v_long_k);
    in.Read("v/diag/v_long_k_0",v_long_k_0);

    // Build k std::vectors
    v_long_k.zeros(path.mag_ks.size());
    for (uint32_t k_i=0; k_i<path.mag_ks.size(); ++k_i) {
      for (uint32_t k_i_v=0; k_i_v<k_v.size(); ++k_i_v) {
        if (fequal(path.mag_ks[k_i],k_v(k_i_v),1.e-8))
          v_long_k(k_i) = tmp_v_long_k(k_i_v);
      }
    }
  }

  // Calculate constants
  uint32_t N1 = path.species_list[species_a_i]->n_part;
  uint32_t N2 = path.species_list[species_b_i]->n_part;
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
double IlkkaPairAction::CalcV(double r, double r_p, const uint32_t level)
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
double IlkkaPairAction::CalcU(double r, double r_p, double s, const uint32_t level)
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
double IlkkaPairAction::CalcdUdBeta(double r, double r_p, double s, const uint32_t level)
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

/// Calculate the V_long value
double IlkkaPairAction::CalcVLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rho_k(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  size_t n_ks = path.ks.size();
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t k_i=0; k_i<n_ks; k_i++) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      tot += v_long_k(k_i)*CMag2(rho_k(path.bead_loop(b_i),species_a_i)(k_i),rho_k(path.bead_loop(b_i),species_b_i)(k_i));
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + v_long_k_0 + v_long_r_0;
}

/// Calculate the ULong value
double IlkkaPairAction::CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
{
  // Get rho k
  field<vec<std::complex<double>>> &rho_k(path.GetRhoK());

  // Sum over k std::vectors
  uint32_t skip = 1<<level;
  double tot = 0.;
  size_t n_ks = path.ks.size();
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t k_i=0; k_i<n_ks; k_i++) {
    for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip) {
      tot += u_long_k(k_i)*CMag2(rho_k(path.bead_loop(b_i),species_a_i)(k_i),rho_k(path.bead_loop(b_i),species_b_i)(k_i));
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
  field<vec<std::complex<double>>> &rho_k(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  size_t n_ks = path.ks.size();
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t k_i=0; k_i<n_ks; k_i++) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      tot += du_long_k(k_i)*CMag2(rho_k(path.bead_loop(b_i),species_a_i)(k_i),rho_k(path.bead_loop(b_i),species_b_i)(k_i));
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + du_long_k_0 + du_long_r_0;
}

vec<double> IlkkaPairAction::CalcGradientU(const uint32_t b_i, const uint32_t b_j, const uint32_t p_i, const uint32_t p_j, const uint32_t level)
{
  // Constants
  vec<double> r, r_p, r_r_p;
  double r_mag, r_p_mag, r_r_p_mag;
  path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,r_mag,r_p_mag,r_r_p_mag,r,r_p,r_r_p);
  double q = 0.5*(r_mag + r_p_mag);
  double x = q + 0.5*r_r_p_mag;
  double y = q - 0.5*r_r_p_mag;

  // Calculate U and (dU/dxp,dU/dyp)
  double u = 0.;
  vec<double> du_dxp_du_dyp(2);
  eval_NUBspline_2d_d_vg(u_xy_spline,x,y,&u,du_dxp_du_dyp.memptr());

  // Do chain rule to get (dU/dx_i,dU/dy_i,dU/dz_i)
  // dU/dx_i = dU/dxp dxp/dx_i + dU/dyp dyp/dx_i
  // and same for y_i and z_i
  //
  // dxp/dx_i = d(q+0.5*s)/dx_i where q = 0.5*(|r| + |r_p|) and s = |r-r_p|
  //                            with r = r_i-r_j and r_p = r_i'-r_j'
  // d|r|/dx_i = d(sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2))/dx_i
  // d|r|/dx_i = (x_i-x_j)/|r|
  // d|r_p|/dx_i = 0
  // d|r-r_p|/dx_i = d(sqrt((x_i-x_j-x_i'+x_j')^2 + (y_i-y_j-y_i'+y_j')^2 + (z_i-z_j-z_i'+z_j')^2))/dx_i
  // d|r-r_p|/dx_i = (x_i-x_j-x_i'+x_j')/|r-r_p|
  // => dxp/dx = 0.5*((x_i-x_j)/|r| + (x_i-x_j-x_i'+x_j')/|r-r_p|)
  //    dyp/dx = 0.5*((x_i-x_j)/|r| - (x_i-x_j-x_i'+x_j')/|r-r_p|)
  // => dU/dx = 0.5*[dU/dxp*((x_i-x_j)/|r| + (x_i-x_j-x_i'+x_j')/|r-r_p|) + dU/dyp*((x_i-x_j)/|r| - (x_i-x_j-x_i'+x_j')/|r-r_p|)]
  // => (dU/dx,dU/dy,dU/dz) = 0.5*[dU/dxp*(r/|r| + (r-r_p)/|r-r_p|) + dU/dyp*(r/|r| - (r-r_p)/|r-r_p|)]

  // Calculate numerical gradient
  vec<double> r_i_r_mag(r/r_mag), r_r_p_i_r_r_p_mag(r_r_p/r_r_p_mag);
  if (r_mag == 0.)
    r_i_r_mag.zeros();
  if (r_r_p_mag == 0.)
    r_r_p_i_r_r_p_mag.zeros();
  vec<double> tot = -0.5*(du_dxp_du_dyp(0)*(r_i_r_mag + r_r_p_i_r_r_p_mag) + du_dxp_du_dyp(1)*(r_i_r_mag - r_r_p_i_r_r_p_mag));

  // Subtract out long range part
  if (use_long_range) {
    // Limits
    SetLimits(r_u_min, r_u_max, r_mag, r_p_mag);

    // Splines
    double tmp_u, tmp_du_dr;
    eval_NUBspline_1d_d_vg(u_long_r_spline,r_mag,&tmp_u,&tmp_du_dr);
    tot -= 0.5*tmp_du_dr*r_i_r_mag;
  }

  return tot;
}

vec<double> IlkkaPairAction::CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
{
  // Should average to 0
  vec<double> tot;
  tot.zeros(path.n_d);
  return tot;
}

vec<double> IlkkaPairAction::CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t p_i, const uint32_t level)
{
  // Get rho k
  field<vec<std::complex<double>>> &rho_k(path.GetRhoK());

  // Sum over k std::vectors
  uint32_t skip = 1<<level;
  vec<double> tot;
  tot.zeros(path.n_d);
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip) {
      vec<std::complex<double>> &rho_k_b(path.GetRhoK(path(species_a_i,p_i,b_i)));
      tot += u_long_k(k_i)*path.ks[k_i]*(rho_k_b(k_i).real()*rho_k(path.bead_loop(b_i),species_b_i)(k_i).imag() - rho_k_b(k_i).imag()*rho_k(path.bead_loop(b_i),species_b_i)(k_i).real());
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot;
}

