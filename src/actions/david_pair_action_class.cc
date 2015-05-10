#include "david_pair_action_class.h"

void DavidPairAction::ReadFile(std::string file_name)
{
  // Load file
  IO in;
  in.Load(file_name);

  // Form std::strings
  std::stringstream tmp_ss;
  tmp_ss << "/u_kj_" << n_order;
  std::string u_kj_str = tmp_ss.str();
  std::stringstream tmp_ss_2;
  tmp_ss_2 << "/du_kj_dbeta_" << n_order;
  std::string du_kj_dbeta_str = tmp_ss_2.str();

  // Read in and create grid
  double r_start, r_end;
  uint32_t n_grid;
  std::string grid_type;
  vec<double> grid_points;
  in.Read(u_kj_str + "/grid/start", r_start);
  in.Read(u_kj_str + "/grid/end", r_end);
  in.Read(u_kj_str + "/grid/n_grid_points", n_grid);
  in.Read(u_kj_str + "/grid/type", grid_type);
  if (grid_type.find("LOG")!=std::string::npos)
    grid = create_log_grid(r_start, r_end, n_grid);
  else if (grid_type.find("LINEAR")!=std::string::npos)
    grid = create_linear_grid(r_start, r_end, n_grid);
  else {
    in.Read(u_kj_str + "/grid/grid_points", grid_points);
    grid = create_general_grid(grid_points.memptr(), n_grid);
  }

  // Read in taus
  n_tau = max_level + 1;
  taus.set_size(n_tau);
  in.Read(u_kj_str + "/taus", taus);
  bool tauFound = 0;
  for (uint32_t tau_i=0; tau_i<n_tau; tau_i++)
    if (abs(taus(tau_i)-path.tau) < 1.0e-6)
      tauFound = 1;
  if (!tauFound) {
    std::cerr << "ERROR: tau of " << path.tau << " not found." << std::endl;
    std::cerr << "Possible taus: " << taus << std::endl;
    exit(1);
  }

  // Read in potential
  vec<double> V(n_grid);
  in.Read("/potential/data", V);

  // Determine number of values for k,j sum
  n_val = 1;
  for (uint32_t i = 1; i <= n_order; ++i)
    n_val += 1+i;

  // Read in u_kj
  cube<double> tmp_u_kj(n_val,n_grid,n_tau);
  in.Read(u_kj_str + "/data", tmp_u_kj);

  // Boundary conditions
  BCtype_d xBC = {NATURAL, NATURAL};

  // Spline u_kj
  cube<double> tmp_u_kj2(n_val+1,n_grid,n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    for (uint32_t grid_i=0; grid_i<n_grid-1; ++grid_i) {
      tmp_u_kj2(0,grid_i,tau_i) = V(grid_i);
      for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
        tmp_u_kj2(val_i,grid_i,tau_i) = tmp_u_kj(val_i-1,grid_i,tau_i);
    }
    //for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
    //  tmp_u_kj2(val_i,n_grid-1,tau_i) = 0.;
  }

  u_kj.set_size(n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    u_kj(tau_i) = create_multi_NUBspline_1d_d(grid, xBC, n_val+1);
    for (uint32_t val_i=0; val_i<n_val+1; ++val_i) {
      vec<double> tmp_v(n_grid);
      for (uint32_t grid_i=0; grid_i<n_grid; ++grid_i)
        tmp_v(grid_i) = tmp_u_kj2(val_i,grid_i,tau_i);
      set_multi_NUBspline_1d_d(u_kj(tau_i), val_i, tmp_v.memptr());
    }
  }

  // Read in du_kj_dbeta
  cube<double> tmpdu_kj_dbeta(n_val,n_grid,n_tau);
  in.Read(du_kj_dbeta_str + "/data", tmpdu_kj_dbeta);

  // Spline du_kj_dbeta
  cube<double> tmpdu_kj_dbeta2(n_val+1,n_grid,n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    for (uint32_t grid_i=0; grid_i<n_grid-1; ++grid_i) {
      tmpdu_kj_dbeta2(0,grid_i,tau_i) = V(grid_i);
      for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
        tmpdu_kj_dbeta2(val_i,grid_i,tau_i) = tmpdu_kj_dbeta(val_i-1,grid_i,tau_i);
    }
    //for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
    //  tmpdu_kj_dbeta2(val_i,n_grid-1,tau_i) = 0.;
  }
  du_kj_dbeta.set_size(n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    du_kj_dbeta(tau_i) = create_multi_NUBspline_1d_d(grid, xBC, n_val+1);
    for (uint32_t val_i=0; val_i<n_val+1; ++val_i) {
      vec<double> tmp_v(n_grid);
      for (uint32_t grid_i=0; grid_i<n_grid; ++grid_i)
        tmp_v(grid_i) = tmpdu_kj_dbeta2(val_i,grid_i,tau_i);
      set_multi_NUBspline_1d_d(du_kj_dbeta(tau_i), val_i, tmp_v.memptr());
    }
  }

  // long range
  if (use_long_range) {
    // Set ks
    double k_v_cut;
    in.Read("long_range/k_cut", k_v_cut);
    uint32_t n_k_v;
    in.Read("long_range/n_k", n_k_v);
    vec<double> k_v(n_k_v);
    in.Read("long_range/k_points", k_v);
    // Read in k space potential
    vec<double> t_v_long_k(n_k_v);
    in.Read("long_range/u_k", t_v_long_k);
    v_long_k = t_v_long_k/path.vol;

    // Build k std::vectors
    u_long_k.zeros(path.mag_ks.size());
    du_long_k.zeros(path.mag_ks.size());
    for (uint32_t k_i_v=0; k_i_v<k_v.size(); ++k_i_v) {
      if (fequal(0.,k_v(k_i_v),1.e-8))
        v_long_k_0 = v_long_k(k_i_v);
      for (uint32_t k_i=0; k_i<path.mag_ks.size(); ++k_i) {
        if (fequal(path.mag_ks[k_i],k_v(k_i_v),1.e-8)) {
          u_long_k(k_i) = v_long_k(k_i_v)*path.tau;
          du_long_k(k_i) = v_long_k(k_i_v);
        }
      }
    }

    // Set constants
    in.Read("squarer/v_image",v_long_r_0);
    u_long_r_0 = v_long_r_0*path.tau; //FIXME: Set to single tau value
    du_long_r_0 = v_long_r_0;
    u_long_k_0 = v_long_k_0*path.tau;
    du_long_k_0 = v_long_k_0;

    // Calculate constants
    uint32_t N1 = path.species_list[species_a_i]->n_part;
    uint32_t N2 = path.species_list[species_b_i]->n_part;
    if (species_a_i == species_b_i) { // homologous
      du_long_k_0 *= 0.5*N1*N1*path.n_bead; //FIXME: Confirm this is correct
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

}

/// Calculate the U(r,r') value when given r and r' and the level 
double DavidPairAction::CalcV(double r, double r_p, const uint32_t level)
{
  // Limits
  double r_min, r_max;
  GetLimits(r_min, r_max, r, r_p, grid);

  // This is the endpoint action
  double v;
  vec<double> r_vals(n_val+1), r_p_vals(n_val+1), q_vals(n_val+1);
  eval_multi_NUBspline_1d_d(u_kj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(u_kj(level),r_p,r_p_vals.memptr());
  v = 0.5*(r_vals(0) + r_p_vals(0));

  return v;
}

/// Calculate the U(r,r') value when given r and r' and the level 
double DavidPairAction::CalcU(double r, double r_p, double s, const uint32_t level)
{
  // Constants
  double q = 0.5*(r + r_p);
  double z = r - r_p;

  // Limits
  double r_min, r_max;
  GetLimits(r_min, r_max, r, r_p, grid);

  // This is the endpoint action
  double u;
  vec<double> r_vals(n_val+1), r_p_vals(n_val+1), q_vals(n_val+1);
  eval_multi_NUBspline_1d_d(u_kj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(u_kj(level),r_p,r_p_vals.memptr());
  u = 0.5*(r_vals(1) + r_p_vals(1));

  // Add in off-diagonal terms
  if (s>0.0 && q<r_max) {
    vec<double> u_q_vals(n_val+1);
    eval_multi_NUBspline_1d_d(u_kj(level),q,u_q_vals.memptr());
    double z_2 = z*z;
    double s_2 = s*s;
    double i_s_2 = 1./s_2;
    double s_2_k = s_2;
    for (uint32_t k=1; k<=n_order; k++) {
      double z_2_j = 1;
      double current_s = s_2_k;
      for (uint32_t j=0; j<=k; j++) {
        // indexing into the 2darray
        double u_cof = u_q_vals(k*(k+1)/2 + (j+1));
        u += (u_cof)*z_2_j*current_s;
        z_2_j *= z_2;
        current_s *= i_s_2;
      }
      s_2_k *= s_2;
    }
  }

  return u;

}

/// Calculate the U(r,r'), dU(r,r'), and V(r,r') value when given r and r' and the level 
double DavidPairAction::CalcdUdBeta(double r, double r_p, double s, const uint32_t level)
{
  // Constants
  double q = 0.5*(r + r_p);
  double z = r - r_p;

  // Limits
  double r_min, r_max;
  GetLimits(r_min, r_max, r, r_p, grid);

  // This is the endpoint action
  double u, v, du;
  vec<double> r_vals(n_val+1), r_p_vals(n_val+1);
  eval_multi_NUBspline_1d_d(u_kj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(u_kj(level),r_p,r_p_vals.memptr());
  v = 0.5*(r_vals(0) + r_p_vals(0));
  u = 0.5*(r_vals(1) + r_p_vals(1));
  eval_multi_NUBspline_1d_d(du_kj_dbeta(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(du_kj_dbeta(level),r_p,r_p_vals.memptr());
  du = 0.5*(r_vals(1) + r_p_vals(1));

  // Compensate for potential, which is subtracted from diaganal action in dm file.
  du += v;

  // Add in off-diagonal terms
  if (s > 0.0 && q<r_max) {
    vec<double> u_q_vals(n_val+1), du_q_vals(n_val+1);
    eval_multi_NUBspline_1d_d(u_kj(level),q,u_q_vals.memptr());
    eval_multi_NUBspline_1d_d(du_kj_dbeta(level),q,du_q_vals.memptr());
    double z_2 = z*z;
    double s_2 = s*s;
    double i_s_2 = 1./s_2;
    double s_2_k = s_2;
    for (uint32_t k=1; k<=n_order; k++) {
      double z_2_j = 1;
      double current_s = s_2_k;
      for (uint32_t j=0; j<=k; j++){
        // indexing into the 2darray
        double u_cof = u_q_vals(k*(k+1)/2+j+1);
        double du_cof = du_q_vals(k*(k+1)/2+j+1);
        u += (u_cof)*z_2_j*current_s;
        du += (du_cof)*z_2_j*current_s;
        z_2_j *= z_2;
        current_s *= i_s_2;
      }
      s_2_k *= s_2;
    }
  }

  return du;
}

/// Calculate the V_long value
double DavidPairAction::CalcVLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += rhok2*v_long_k(k_i);
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + v_long_k_0 + v_long_r_0;
}

/// Calculate the ULong value
double DavidPairAction::CalcULong(const uint32_t b_0, const uint32_t b_1, const uint32_t level)
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  uint32_t skip = 1<<level;
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint32_t b_i=b_0; b_i<b_1; b_i+=skip) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += u_long_k(k_i)*rhok2;
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot;
}

/// Calculate the dUdBetaLong value
double DavidPairAction::CalcdUdBetaLong()
{
  // Get rho k
  field<vec<std::complex<double>>> &rhoK(path.GetRhoK());

  // Sum over k std::vectors
  double tot = 0.;
  for (uint32_t k_i=0; k_i<path.ks.size(); k_i++) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
      double rhok2 = CMag2(rhoK(path.bead_loop(b_i),species_a_i)(k_i),rhoK(path.bead_loop(b_i),species_b_i)(k_i));
      tot += du_long_k(k_i)*rhok2;
    }
  }

  if (species_b_i != species_a_i)
    tot *= 2.;

  return tot + du_long_k_0 + du_long_r_0;
}
