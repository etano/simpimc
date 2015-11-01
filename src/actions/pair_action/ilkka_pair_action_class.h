#ifndef SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Pair action class derived from Ilkka Kylanpaa's squarer
class IlkkaPairAction : public PairAction
{
private:
  double du_long_k_0; ///< Constant correction to the k components of the beta derivative of the long range pair action
  double du_long_r_0; ///< Constant correction to spatial beta derivative of the long range pair action
  double r_u_max; ///< Maximum distance possible in the pair action spline
  double r_u_min; ///< Minimum distance possible in the pair action spline
  double r_v_max; ///< Maximum distance possible in the pair potential spline
  double r_v_min; ///< Minimum distance possible in the pair potential spline
  double r_du_max; ///< Maximum distance possible in the beta derivative of the pair action spline
  double r_du_min; ///< Minimum distance possible in the beta derivative of the pair action spline
  double r_v_long_max; ///< Maximum distance possible in the long range pair potential spline
  double r_v_long_min; ///< Minimum distance possible in the long range pair potential spline
  double u_long_k_0; ///< Constant correction to the k components of the long range pair action
  double u_long_r_0; ///< Constant correction to spatial long range pair action
  double v_long_k_0; ///< Constant correction to the k components of the long range pair potential
  double v_long_r_0; ///< Constant correction to spatial long range pair potential
  vec<double> u_long_k; ///< Vector of k components of the long range pair action
  vec<double> du_long_k; ///< Vector of k components of the beta derivative of the long range pair action
  vec<double> v_long_k; ///< Vector of k components of the long range pair potential
  NUBspline_1d_d *u_long_r_spline; ///< Spline of the long range pair action
  NUBspline_1d_d *du_long_r_spline; ///< Spline of the beta derivative of the long range pair action
  NUBspline_1d_d *v_r_spline; ///< Spline of the pair potential
  NUBspline_1d_d *v_long_r_spline; ///< Spline of the long range pair potential
  NUBspline_2d_d *u_xy_spline; ///< Spline of the pair action
  NUBspline_2d_d *du_xy_spline; ///< Spline of the beta derivative of the pair action

  /// Calculate the potential
  virtual double CalcV(double r, double r_p, const uint32_t level)
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

  /// Calculate the long ranged part of the potential
  virtual double CalcVLong()
  {
    // Get rho k
    field<vec<std::complex<double>>> &rho_k_a(species_a->GetRhoK());
    field<vec<std::complex<double>>> &rho_k_b(species_b->GetRhoK());

    // Sum over k std::vectors
    double tot = 0.;
    size_t n_ks = path.ks.vecs.size();
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t k_i=0; k_i<n_ks; k_i++)
      for (uint32_t b_i=0; b_i<path.GetNBead(); b_i++)
        tot += v_long_k(k_i)*CMag2(rho_k_a(b_i)(k_i),rho_k_b(b_i)(k_i));

    if (species_a != species_b)
      tot *= 2.;

    return tot + v_long_k_0 + v_long_r_0;
  }

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level)
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

  /// Calculate the long ranged part of the action
  virtual double CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
  {
    // Get rho k
    field<vec<std::complex<double>>> &rho_k_a(species_a->GetRhoK());
    field<vec<std::complex<double>>> &rho_k_b(species_b->GetRhoK());

    // Sum over k std::vectors
    uint32_t skip = 1<<level;
    double tot = 0.;
    size_t n_ks = path.ks.vecs.size();
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t k_i=0; k_i<n_ks; k_i++)
      for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip)
        tot += u_long_k(k_i)*CMag2(rho_k_a(species_a->bead_loop(b_i))(k_i),rho_k_b(species_b->bead_loop(b_i))(k_i));

    if (species_a != species_b)
      tot *= 2.;

    return tot;
  }

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level)
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

  /// Calculate the long ranged part of the beta derivative of the action
  virtual double CalcdUdBetaLong()
  {
    // Get rho k
    field<vec<std::complex<double>>> &rho_k_a(species_a->GetRhoK());
    field<vec<std::complex<double>>> &rho_k_b(species_b->GetRhoK());

    // Sum over k std::vectors
    double tot = 0.;
    size_t n_ks = path.ks.vecs.size();
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t k_i=0; k_i<n_ks; k_i++)
      for (uint32_t b_i=0; b_i<path.GetNBead(); b_i++)
        tot += du_long_k(k_i)*CMag2(rho_k_a(species_a->bead_loop(b_i))(k_i),rho_k_b(species_b->bead_loop(b_i))(k_i));

    if (species_a != species_b)
      tot *= 2.;

    return tot + du_long_k_0 + du_long_r_0;
  }

  /// Calculate the gradient of the action for the particle pair p_i, p_j in the direction of particle p_i
  virtual vec<double> CalcGradientU(const uint32_t b_i, const uint32_t b_j, const uint32_t p_i, const uint32_t p_j, const uint32_t level)
  {
    // Constants
    vec<double> r, r_p, r_r_p;
    double r_mag, r_p_mag, r_r_p_mag;
    path.DrDrpDrrp(b_i,b_j,species_a,species_b,p_i,p_j,r_mag,r_p_mag,r_r_p_mag,r,r_p,r_r_p);
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
    //           = (x_i-x_j)/|r|
    // d|r_p|/dx_i = 0
    // d|r-r_p|/dx_i = d(sqrt((x_i-x_j-x_i'+x_j')^2 + (y_i-y_j-y_i'+y_j')^2 + (z_i-z_j-z_i'+z_j')^2))/dx_i
    //               = (x_i-x_j-x_i'+x_j')/|r-r_p|
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

  /// Calculate the gradient of the long ranged part of the action for all particles
  virtual vec<double> CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
  {
    // Should average to 0
    return zeros<vec<double>>(path.GetND());
  }

  /// Calculate the gradient of the long ranged part of the action in the direction of p_i
  virtual vec<double> CalcGradientULong(const uint32_t b_0, const uint32_t b_1, const uint32_t p_i, const uint32_t level)
  {
    // Get rho k
    field<vec<std::complex<double>>> &rho_k_b(species_b->GetRhoK());

    // Sum over k std::vectors
    uint32_t skip = 1<<level;
    vec<double> tot(zeros<vec<double>>(path.GetND()));
    for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip) {
      vec<std::complex<double>> &rho_k_a_i(species_a->GetBead(p_i,b_i)->GetRhoK());
      for (uint32_t k_i=0; k_i<path.ks.vecs.size(); k_i++) {
        tot += u_long_k(k_i)*path.ks.vecs[k_i]*(rho_k_a_i(k_i).real()*rho_k_b(species_b->bead_loop(b_i))(k_i).imag() - rho_k_a_i(k_i).imag()*rho_k_b(species_b->bead_loop(b_i))(k_i).real());
      }
    }

    if (species_a != species_b)
      tot *= 2.;

    return tot;
  }

public:
  /// Constructor only calls Init
  IlkkaPairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    // Get file
    std::string file_name = in.GetAttribute<std::string>("file");
    out.Write("Actions/"+name+"/file", file_name);

    // Load file
    IO pa_in;
    pa_in.Load(file_name);

    // Boundary conditions
    BCtype_d xBC = {NATURAL, NATURAL}; // HACK: Is this correct?

    // Read in u
    uint32_t n_x_u, n_y_u;
    pa_in.Read("u/off_diag/n_x", n_x_u);
    pa_in.Read("u/off_diag/n_y", n_y_u);
    vec<double> x_u(n_x_u);
    vec<double> y_u(n_y_u);
    mat<double> u_xy(n_x_u,n_y_u);
    pa_in.Read("u/off_diag/x", x_u);
    pa_in.Read("u/off_diag/y", y_u);
    pa_in.Read("u/off_diag/u_xy", u_xy);

    // Spline u
    NUgrid* x_u_grid = create_general_grid(x_u.memptr(), x_u.size());
    NUgrid* y_u_grid = create_general_grid(y_u.memptr(), y_u.size());
    u_xy_spline = create_NUBspline_2d_d(x_u_grid, y_u_grid, xBC, xBC, u_xy.memptr());

    // u long range
    if (use_long_range) {
      // Read in r
      uint32_t n_r_u;
      pa_in.Read("u/diag/n_r_long", n_r_u);
      vec<double> r_u(n_r_u);
      vec<double> u_long_r(n_r_u);
      pa_in.Read("u/diag/r_long", r_u);
      pa_in.Read("u/diag/u_long_r", u_long_r);
      pa_in.Read("u/diag/u_long_r_0",u_long_r_0);

      // Spline r
      NUgrid* r_u_grid = create_general_grid(r_u.memptr(), r_u.size());
      r_u_min = r_u_grid->start;
      r_u_max = r_u_grid->end;
      u_long_r_spline = create_NUBspline_1d_d(r_u_grid, xBC, u_long_r.memptr());

      // Read in k
      uint32_t n_k_u;
      pa_in.Read("u/diag/n_k", n_k_u);
      vec<double> k_u(n_k_u);
      vec<double> tmp_u_long_k(n_k_u);
      pa_in.Read("u/diag/k",k_u);
      pa_in.Read("u/diag/u_long_k",tmp_u_long_k);
      pa_in.Read("u/diag/u_long_k_0",u_long_k_0);

      // Build k std::vectors
      u_long_k.zeros(path.ks.mags.size());
      for (uint32_t k_i=0; k_i<path.ks.mags.size(); ++k_i) {
        for (uint32_t k_i_u=0; k_i_u<k_u.size(); ++k_i_u) {
          if (fequal(path.ks.mags[k_i],k_u(k_i_u),1.e-8))
            u_long_k(k_i) = tmp_u_long_k(k_i_u);
        }
      }
    }

    // Read in du
    uint32_t n_x_du, n_y_du;
    pa_in.Read("du/off_diag/n_x", n_x_du);
    pa_in.Read("du/off_diag/n_y", n_y_du);
    vec<double> x_du(n_x_du);
    vec<double> y_du(n_y_du);
    mat<double> du_xy(n_x_du,n_y_du);
    pa_in.Read("du/off_diag/x", x_du);
    pa_in.Read("du/off_diag/y", y_du);
    pa_in.Read("du/off_diag/du_xy", du_xy);

    // Spline du
    NUgrid* x_du_grid = create_general_grid(x_du.memptr(), x_du.size());
    NUgrid* y_du_grid = create_general_grid(y_du.memptr(), y_du.size());
    du_xy_spline = create_NUBspline_2d_d(x_du_grid, y_du_grid, xBC, xBC, du_xy.memptr());

    // du long range
    if (use_long_range) {
      // Read in r
      uint32_t n_r_du;
      pa_in.Read("du/diag/n_r_long", n_r_du);
      vec<double> r_du(n_r_du);
      vec<double> du_long_r(n_r_du);
      pa_in.Read("du/diag/r_long", r_du);
      pa_in.Read("du/diag/du_long_r", du_long_r);
      pa_in.Read("du/diag/du_long_r_0",du_long_r_0);

      // Spline r
      NUgrid* r_du_grid = create_general_grid(r_du.memptr(), r_du.size());
      r_du_min = r_du_grid->start;
      r_du_max = r_du_grid->end;
      du_long_r_spline = create_NUBspline_1d_d(r_du_grid, xBC, du_long_r.memptr());

      // Read in k
      uint32_t n_k_du;
      pa_in.Read("du/diag/n_k", n_k_du);
      vec<double> k_du(n_k_du);
      vec<double> tmp_du_long_k(n_k_du);
      pa_in.Read("du/diag/k",k_du);
      pa_in.Read("du/diag/du_long_k",tmp_du_long_k);
      pa_in.Read("du/diag/du_long_k_0",du_long_k_0);

      // Build k std::vectors
      du_long_k.zeros(path.ks.mags.size());
      for (uint32_t k_i=0; k_i<path.ks.mags.size(); ++k_i) {
        for (uint32_t k_i_du=0; k_i_du<k_du.size(); ++k_i_du) {
          if (fequal(path.ks.mags[k_i],k_du(k_i_du),1.e-8))
            du_long_k(k_i) = tmp_du_long_k(k_i_du);
        }
      }

    }

    // Read in v
    uint32_t n_r_v;
    pa_in.Read("v/diag/n_r", n_r_v);
    vec<double> r_v(n_r_v);
    vec<double> v_r(n_r_v);
    pa_in.Read("v/diag/r", r_v);
    pa_in.Read("v/diag/v_r", v_r);

    // Spline v
    NUgrid* r_v_grid = create_general_grid(r_v.memptr(), r_v.size());
    r_v_min = r_v_grid->start;
    r_v_max = r_v_grid->end;
    v_r_spline = create_NUBspline_1d_d(r_v_grid, xBC, v_r.memptr());

    // v long range
    if (use_long_range) {
      // Read in r
      uint32_t n_r_v_long;
      pa_in.Read("v/diag/n_r_long", n_r_v_long);
      vec<double> r_v_long(n_r_v_long);
      vec<double> v_long_r(n_r_v_long);
      pa_in.Read("v/diag/r_long",r_v_long);
      pa_in.Read("v/diag/v_long_r",v_long_r);
      pa_in.Read("v/diag/v_long_r_0",v_long_r_0);

      // Spline r
      NUgrid* r_v_long_grid = create_general_grid(r_v_long.memptr(), r_v_long.size());
      r_v_long_min = r_v_long_grid->start;
      r_v_long_max = r_v_long_grid->end;
      v_long_r_spline = create_NUBspline_1d_d(r_v_long_grid, xBC, v_long_r.memptr());

      // Read in k
      uint32_t n_k_v;
      pa_in.Read("v/diag/n_k", n_k_v);
      vec<double> k_v(n_k_v);
      vec<double> tmp_v_long_k(n_k_v);
      pa_in.Read("v/diag/k",k_v);
      pa_in.Read("v/diag/v_long_k",tmp_v_long_k);
      pa_in.Read("v/diag/v_long_k_0",v_long_k_0);

      // Build k std::vectors
      v_long_k.zeros(path.ks.mags.size());
      for (uint32_t k_i=0; k_i<path.ks.mags.size(); ++k_i) {
        for (uint32_t k_i_v=0; k_i_v<k_v.size(); ++k_i_v) {
          if (fequal(path.ks.mags[k_i],k_v(k_i_v),1.e-8))
            v_long_k(k_i) = tmp_v_long_k(k_i_v);
        }
      }
    }

    // Calculate constants
    if (species_a == species_b) { // homologous
      du_long_k_0 *= 0.5*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      du_long_r_0 *= -0.5*species_a->GetNPart()*path.GetNBead();
      v_long_k_0 *= 0.5*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      v_long_r_0 *= -0.5*species_a->GetNPart()*path.GetNBead();
    } else { // heterologous
      du_long_k_0 *= species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      du_long_r_0 *= 0.;
      v_long_k_0 *= species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      v_long_r_0 *= 0.;
    }

  }

};

#endif // SIMPIMC_ACTIONS_ILKKA_PAIR_ACTION_CLASS_H_
