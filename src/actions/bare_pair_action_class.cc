#include "bare_pair_action_class.h"

void BarePairAction::ReadFile(std::string file_name)
{
  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Load file
  IO in;
  in.Load(file_name);

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
    v_long_k_0 *= 0.5*N1*N1*path.n_bead;
    v_long_r_0 *= -0.5*N1*path.n_bead;
  } else { // heterologous
    v_long_k_0 *= N1*N2*path.n_bead;
    v_long_r_0 *= 0.*path.n_bead;
  }

}

/// Calculate the V(r,r') value when given r and r' and the level 
double BarePairAction::CalcV(double r, double r_p, const uint32_t level)
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
double BarePairAction::CalcU(double r, double r_p, double s, const uint32_t level)
{
  uint32_t skip = 1>>level;
  double level_tau = skip*path.tau;
  return level_tau*CalcV(r,r_p,level);
}

/// Calculate the dU(r,r') value when given r and r' and the level 
double BarePairAction::CalcdUdBeta(double r, double r_p, double s, const uint32_t level)
{
  return CalcV(r,r_p,level);
}

/// Calculate the ULong value
double BarePairAction::CalcVLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    if (path.mag_ks[k_i] < k_cut) {
      for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
        double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
        tot += rhok2*v_long_k(k_i);
      }
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + v_long_k_0 + v_long_r_0;
}

/// Calculate the ULong value
double BarePairAction::CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  uint32_t skip = 1<<level;
  double tot = 0.;
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    if (path.mag_ks[k_i] < k_cut) {
      for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip) {
        double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
        tot += v_long_k(k_i)*rhok2;
      }
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
double BarePairAction::CalcdUdBetaLong()
{
  return CalcVLong();
}
