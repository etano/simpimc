#include "david_pair_action_class.h"

void DavidPairAction::ReadFile(std::string file_name)
{
  // Load file
  IO in;
  in.Load(file_name);

  // Form std::strings
  std::stringstream tmp_ss;
  tmp_ss << "/ukj" << n_order;
  std::string ukj_str = tmp_ss.str();
  std::stringstream tmp_ss_2;
  tmp_ss_2 << "/dukj_dbeta" << n_order;
  std::string dukj_dbeta_str = tmp_ss_2.str();

  // Read in and create grid
  double r_start, r_end;
  uint32_t n_grid;
  std::string grid_type;
  vec<double> grid_points;
  in.Read(ukj_str + "/Grid/Start", r_start);
  in.Read(ukj_str + "/Grid/End", r_end);
  in.Read(ukj_str + "/Grid/NumGridPoints", n_grid);
  in.Read(ukj_str + "/Grid/Type", grid_type);
  if (grid_type == "LOG")
    grid = create_log_grid(r_start, r_end, n_grid);
  else if (grid_type == "LINEAR")
    grid = create_linear_grid(r_start, r_end, n_grid);
  else {
    in.Read(ukj_str + "/Grid/GridPoints", grid_points);
    grid = create_general_grid(grid_points.memptr(), n_grid);
  }

  // Read in taus
  n_tau = max_level + 1;
  taus.set_size(n_tau);
  in.Read(ukj_str + "/Taus", taus);
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
  in.Read("/Potential/Data", V);

  // Determine number of values for k,j sum
  n_val = 1;
  for (uint32_t i = 1; i <= n_order; ++i)
    n_val += 1+i;

  // Read in ukj
  cube<double> tmp_ukj(n_val,n_grid,n_tau);
  in.Read(ukj_str + "/Data", tmp_ukj);

  // Boundary conditions
  BCtype_d xBC = {NATURAL, FLAT}; // HACK: Is this correct?

  // Spline ukj
  cube<double> tmp_ukj2(n_val+1,n_grid,n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    for (uint32_t grid_i=0; grid_i<n_grid-1; ++grid_i) {
      tmp_ukj2(0,grid_i,tau_i) = V(grid_i);
      for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
        tmp_ukj2(val_i,grid_i,tau_i) = tmp_ukj(val_i-1,grid_i,tau_i);
    }
    //for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
    //  tmp_ukj2(val_i,n_grid-1,tau_i) = 0.;
  }

  ukj.set_size(n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    ukj(tau_i) = create_multi_NUBspline_1d_d(grid, xBC, n_val+1);
    for (uint32_t val_i=0; val_i<n_val+1; ++val_i) {
      vec<double> tmp_v(n_grid);
      for (uint32_t grid_i=0; grid_i<n_grid; ++grid_i)
        tmp_v(grid_i) = tmp_ukj2(val_i,grid_i,tau_i);
      set_multi_NUBspline_1d_d(ukj(tau_i), val_i, tmp_v.memptr());
    }
  }

  // Read in dukj_dbeta
  cube<double> tmpdukj_dbeta(n_val,n_grid,n_tau);
  in.Read(dukj_dbeta_str + "/Data", tmpdukj_dbeta);

  // Spline dukj_dbeta
  cube<double> tmpdukj_dbeta2(n_val+1,n_grid,n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    for (uint32_t grid_i=0; grid_i<n_grid-1; ++grid_i) {
      tmpdukj_dbeta2(0,grid_i,tau_i) = V(grid_i);
      for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
        tmpdukj_dbeta2(val_i,grid_i,tau_i) = tmpdukj_dbeta(val_i-1,grid_i,tau_i);
    }
    //for (uint32_t val_i=1; val_i<n_val+1; ++val_i)
    //  tmpdukj_dbeta2(val_i,n_grid-1,tau_i) = 0.;
  }
  dukj_dbeta.set_size(n_tau);
  for(uint32_t tau_i=0; tau_i<n_tau; tau_i++) {
    dukj_dbeta(tau_i) = create_multi_NUBspline_1d_d(grid, xBC, n_val+1);
    for (uint32_t val_i=0; val_i<n_val+1; ++val_i) {
      vec<double> tmp_v(n_grid);
      for (uint32_t grid_i=0; grid_i<n_grid; ++grid_i)
        tmp_v(grid_i) = tmpdukj_dbeta2(val_i,grid_i,tau_i);
      set_multi_NUBspline_1d_d(dukj_dbeta(tau_i), val_i, tmp_v.memptr());
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
  vec<double> r_vals(n_val+1), r_p_vals(n_val+1), qVals(n_val+1);
  eval_multi_NUBspline_1d_d(ukj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(ukj(level),r_p,r_p_vals.memptr());
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
  vec<double> r_vals(n_val+1), r_p_vals(n_val+1), qVals(n_val+1);
  eval_multi_NUBspline_1d_d(ukj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(ukj(level),r_p,r_p_vals.memptr());
  u = 0.5*(r_vals(1) + r_p_vals(1));

  // Add in off-diagonal terms
  if (s>0.0 && q<r_max) {
    vec<double> UqVals(n_val+1);
    eval_multi_NUBspline_1d_d(ukj(level),q,UqVals.memptr());
    double z2 = z*z;
    double s2 = s*s;
    double s2inverse = 1./s2;
    double Sto2k = s2;
    for (uint32_t k=1; k<=n_order; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (uint32_t j=0; j<=k; j++) {
        // indexing into the 2darray
        double Ucof = UqVals(k*(k+1)/2 + (j+1));
        u += (Ucof)*Zto2j*currS;
        Zto2j *= z2;
        currS *= s2inverse;
      }
      Sto2k *= s2;
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
  eval_multi_NUBspline_1d_d(ukj(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(ukj(level),r_p,r_p_vals.memptr());
  v = 0.5*(r_vals(0) + r_p_vals(0));
  u = 0.5*(r_vals(1) + r_p_vals(1));
  eval_multi_NUBspline_1d_d(dukj_dbeta(level),r,r_vals.memptr());
  eval_multi_NUBspline_1d_d(dukj_dbeta(level),r_p,r_p_vals.memptr());
  du = 0.5*(r_vals(1) + r_p_vals(1));

  // Compensate for potential, which is subtracted from diaganal action in dm file.
  du += v;

  // Add in off-diagonal terms
  if (s > 0.0 && q<r_max) {
    vec<double> UqVals(n_val+1), dUqVals(n_val+1);
    eval_multi_NUBspline_1d_d(ukj(level),q,UqVals.memptr());
    eval_multi_NUBspline_1d_d(dukj_dbeta(level),q,dUqVals.memptr());
    double z2 = z*z;
    double s2 = s*s;
    double s2inverse = 1./s2;
    double Sto2k = s2;
    for (uint32_t k=1; k<=n_order; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (uint32_t j=0; j<=k; j++){
        // indexing into the 2darray
        double Ucof = UqVals(k*(k+1)/2+j+1);
        double dUcof = dUqVals(k*(k+1)/2+j+1);
        u += (Ucof)*Zto2j*currS;
        du += (dUcof)*Zto2j*currS;
        Zto2j *= z2;
        currS *= s2inverse;
      }
      Sto2k *= s2;
    }
  }

  return du;
}
