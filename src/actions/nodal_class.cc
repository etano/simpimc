#include "nodal_class.h"

// Initialize parameters
void Nodal::Init(Input &in)
{
  // Read in things
  n_images = in.GetAttribute<int>("n_images",0);
  is_importance_weight = in.GetAttribute<bool>("is_importance_weight",false);

  // Set species variables
  species = in.GetAttribute<std::string>("species");
  species_list.push_back(species);
  std::cout << "Setting up nodal action for " << species << "..." << std::endl;
  max_level = in.GetAttribute<uint32_t>("max_level",0);
  path.GetSpeciesInfo(species,species_i);
  n_part = path.species_list[species_i]->n_part;
  i_4_lambda_tau = 1./(4.*path.species_list[species_i]->lambda*path.tau);
  path.species_list[species_i]->fixed_node = true;

  // Nodal distance things
  use_nodal_distance = in.GetAttribute<bool>("use_nodal_distance",false);
  if (use_nodal_distance) {
    std::string dist_type_name = in.GetAttribute<std::string>("dist_type");
    std::cout << "Setting up " << dist_type_name << " distance measure for nodal action..." << std::endl;
    if (dist_type_name == "NewtonRaphson") {
      dist_type = 0;
      max_dist_steps = in.GetAttribute<int>("max_dist_steps",1);
    } else if (dist_type_name == "LineSearch")
      dist_type = 1;
    else if (dist_type_name == "Max")
      dist_type = 2;
    else if (dist_type_name == "Hybrid")
      dist_type = 3;
    else if (dist_type_name == "Hyperplane")
      dist_type = 4;
    else {
      std::cerr << "ERROR: Unknown dist_type!" << std::endl;
      exit(1);
    }
    dist_tolerance = in.GetAttribute<double>("tolerance",1.e-4);
    dist.zeros(path.n_bead);
    dist_c.zeros(path.n_bead);
  }

  // Write things to file
  out.Write("Actions/"+name+"/n_images", n_images);
  out.Write("Actions/"+name+"/species", species);
  out.Write("Actions/"+name+"/max_level", max_level);

  // Setup splines
  SetupSpline();

  // Set up determinants
  rho_f.zeros(path.n_bead);
  rho_f_c.zeros(path.n_bead);
  grad_rho_f.set_size(path.n_bead,n_part);
  grad_rho_f_c.set_size(path.n_bead,n_part);
  for (uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
    for (uint32_t p_i=0; p_i<n_part; ++p_i) {
      grad_rho_f(b_i,p_i).zeros(path.n_d);
      grad_rho_f_c(b_i,p_i).zeros(path.n_d);
    }
  }

  // Test initial configuration
  is_first_time.resize(path.n_bead);
  for (uint32_t b_i=0; b_i<path.n_bead; ++b_i)
    is_first_time[b_i] = true;
  std::vector< std::pair<uint32_t,uint32_t>> particles;
  particles.push_back(std::make_pair(species_i,0));
  bool init_good = 1;
  ModeType t_mode(path.GetMode());
  path.SetMode(NEW_MODE);
  if (GetAction(0, path.n_bead, particles, 0) >= 1.e100) {
    std::cout << "Warning: initializing with broken nodes!" << std::endl;
    init_good = 0;
  }
  start_b = 0;
  end_b = path.n_bead;
  Accept();
  path.SetMode(t_mode);
  out.Write("Actions/"+name+"/init_good", init_good);
}

double Nodal::DActionDBeta()
{
  // Zero if not using the nodal distance or using it as an importance weight
  if (!use_nodal_distance)
    return 0.;

  // Initialize other beads
  std::vector<std::shared_ptr<Bead>> ref_b(n_part);
  std::vector<std::vector<std::shared_ptr<Bead>>> other_b(path.n_bead);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    ref_b[p_i] = path(species_i,p_i,path.ref_bead);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    other_b[0].push_back(path.GetPrevBead(ref_b[p_i],path.ref_bead));

  std::vector<uint32_t> b_i_vec;
  b_i_vec.push_back(other_b[0][0]->b);
  for (uint32_t b_i=1; b_i<path.n_bead; ++b_i) {
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_b[b_i].push_back(path.GetNextBead(other_b[b_i-1][p_i],1));
    b_i_vec.push_back(other_b[b_i][0]->b);
  }

  // Decide which type of action to compute
  return DDistanceActionDBeta(b_i_vec, ref_b, other_b);
}

double Nodal::DDistanceActionDBeta(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b) // FIXME: No need for this extra function
{
  // Compute action from nodal distance
  double tot = 0.;
  double i_lambda_tau = 1./(path.species_list[species_i]->lambda*path.tau);
  for (uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
    // Set slice variables
    uint32_t b_j = path.bead_loop(b_i + 1);
    bool b_i_is_ref(b_i == path.ref_bead);
    bool b_j_is_ref(b_j == path.ref_bead);

    // Compute dist2
    double dist_b_i, dist_b_j, dist2;
    if (!b_i_is_ref)
      dist_b_i = dist(b_i); //GetNodalDistance(b_i_vec[b_i],ref_b,other_b[b_i]);
    if (!b_j_is_ref)
      dist_b_j = dist(b_j); //GetNodalDistance(b_i_vec[b_j],ref_b,other_b[b_j]);
    if (b_i_is_ref || (dist_b_i==0.))
      dist2 = dist_b_j*dist_b_j;
    else if (b_j_is_ref || (dist_b_j==0.))
      dist2 = dist_b_i*dist_b_i;
    else
      dist2 = dist_b_i*dist_b_j;

    // Compute action
    double dist2_i_lambda_tau = dist2*i_lambda_tau;
    double exp_dist2_i_lambda_tau_m1 = expm1(dist2_i_lambda_tau);
    if (std::fpclassify(exp_dist2_i_lambda_tau_m1) == FP_NORMAL)
      tot += dist2_i_lambda_tau/(path.tau*exp_dist2_i_lambda_tau_m1);
  }

  return tot/path.n_bead;
}

double Nodal::GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  // No need to check old action if not using the nodal distance
  if (path.GetMode()==OLD_MODE && !use_nodal_distance)
    return 0.;

  // Decide whether or not to check the node
  if (level > max_level || !path.species_list[species_i]->fermi)
    return 0.;
  bool check_node = false;
  for (auto& p: particles) {
    if (p.first == species_i) {
      check_node = true;
      break;
    }
  }
  if (!check_node)
    return 0.;

  // See if ref slice included
  bool check_all = false;
  if (b1 < path.n_bead)
    check_all = ((b0 <= path.ref_bead) && (b1 >= path.ref_bead));
  else
    check_all = (path.bead_loop(b1) >= path.ref_bead);

  // Set start and end beads
  if (check_all) {
    start_b = 0;
    end_b = path.n_bead-1;
  } else {
    start_b = b0;
    end_b = b1;
  }

  // Set up bead index vector
  uint32_t skip = 1<<level;
  int n_bead_in_move = 1 + (end_b - start_b)/skip;
  std::vector<uint32_t> b_i_vec;
  for (uint32_t b_i=start_b; b_i<=end_b; b_i+=skip)
    b_i_vec.push_back(path.bead_loop(b_i));

  // Initialize other beads
  std::vector<std::shared_ptr<Bead>> ref_b(n_part);
  std::vector<std::vector<std::shared_ptr<Bead>>> other_b(n_bead_in_move);
  int slice_diff_0 = path.bead_loop(start_b) - path.ref_bead;
  int abs_slice_diff_0 = abs(slice_diff_0);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    ref_b[p_i] = path(species_i,p_i,path.ref_bead);
  if (slice_diff_0 >= 0) {
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_b[0].push_back(path.GetNextBead(ref_b[p_i],abs_slice_diff_0)); // FIXME: Perhaps this could be faster
  } else {
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_b[0].push_back(path.GetPrevBead(ref_b[p_i],abs_slice_diff_0));
  }
  for (uint32_t b_i=1; b_i<n_bead_in_move; ++b_i)
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_b[b_i].push_back(path.GetNextBead(other_b[b_i-1][p_i],skip));

  // Decide which type of action to compute
  if (use_nodal_distance)
    return DistanceAction(b_i_vec, ref_b, other_b, n_bead_in_move, check_all);
  else
    return SimpleAction(b_i_vec, ref_b, other_b, n_bead_in_move, check_all);
}

double Nodal::SimpleAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all)
{
  // Compute action
  double tot = 0.;
  if (check_all) {
    std::atomic_bool abort(false);
    #pragma omp parallel for reduction(+:tot) schedule(dynamic) shared(abort) // FIXME: Could be optimized probably
    for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
      if (!abort && b_i_vec[b_i] != path.ref_bead) {
        SetRhoF(b_i_vec[b_i],ref_b,other_b[b_i]);
        if (rho_f(b_i_vec[b_i]) < 0.) {
          tot += 1.e100;
          abort = true;
        }
      }
    }
  } else {
    for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
      if (b_i_vec[b_i] != path.ref_bead)  {
        SetRhoF(b_i_vec[b_i],ref_b,other_b[b_i]);
        if (rho_f(b_i_vec[b_i]) < 0.) {
          tot += 1.e100;
          break;
        }
      }
    }
  }

  return tot;
}

// Form rho_f
void Nodal::SetRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  // Set slice distance
  uint32_t abs_slice_diff = abs(path.bead_loop(b_i)-path.ref_bead);
  uint32_t min_slice_diff = std::min(abs_slice_diff, path.n_bead-abs_slice_diff);

  // Compute Gij matrix
  mat<double> g_ij(n_part,n_part);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    for (uint32_t p_j=0; p_j<n_part; ++p_j)
      g_ij(p_i,p_j) = GetGij(path.Dr(ref_b[p_i], other_b_i[p_j]), min_slice_diff);

  // Compute determinant
  rho_f(b_i) = det(g_ij);

  return;
}

double Nodal::DistanceAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all)
{
  // Set nodal distances
  double tot = 0.;
  vec<double> t_dist(n_bead_in_move);
  if (check_all) {
    std::atomic_bool abort(false);
    #pragma omp parallel for reduction(+:tot) schedule(dynamic) shared(abort) // TODO: Could be optimized probably
    for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
      if (!abort && b_i_vec[b_i] != path.ref_bead)  {
        t_dist(b_i) = GetNodalDistance(b_i_vec[b_i],ref_b,other_b[b_i]);
        if (t_dist(b_i) < 0.) {
          tot += 1.e100;
          abort = true;
        }
      }
    }
  } else {
    for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
      if (b_i_vec[b_i] != path.ref_bead)  {
        t_dist(b_i) = GetNodalDistance(b_i_vec[b_i],ref_b,other_b[b_i]);
        if (t_dist(b_i) < 0.) {
          tot += 1.e100;
          break;
        }
      }
    }
  }

  // Compute action from nodal distance
  if (tot == 0.) {
    uint32_t skip = b_i_vec[1] - b_i_vec[0];
    double i_lambda_level_tau = 1./(path.species_list[species_i]->lambda*skip*path.tau);
    for (uint32_t b_i=0; b_i<n_bead_in_move-1; ++b_i) {

      // Set slice variables
      uint32_t b_j = b_i + 1;
      bool b_i_is_ref = b_i_vec[b_i] == path.ref_bead;
      bool b_j_is_ref = b_i_vec[b_j] == path.ref_bead;

      // Compute action
      if (!b_i_is_ref && (t_dist(b_i)<0.)) {
        tot = 1.e100;
        break;
      } else if (!b_j_is_ref && (t_dist(b_j)<0.)) {
        tot = 1.e100;
        break;
      } else if (b_i_is_ref || (t_dist(b_i)==0.))
        tot -= log1p(-path.FastExp(-t_dist(b_j)*t_dist(b_j)*i_lambda_level_tau));
      else if (b_j_is_ref || (t_dist(b_j)==0.))
        tot -= log1p(-path.FastExp(-t_dist(b_i)*t_dist(b_i)*i_lambda_level_tau));
      else
        tot -= log1p(-path.FastExp(-t_dist(b_i)*t_dist(b_j)*i_lambda_level_tau));

      // Store if first time or in new mode
      if (path.GetMode()==NEW_MODE || is_first_time[b_i_vec[b_i]]) {
        dist(b_i_vec[b_i]) = t_dist(b_i);
        is_first_time[b_i_vec[b_i]] = false;
      }
      if (path.GetMode()==NEW_MODE || is_first_time[b_i_vec[b_j]]) {
        dist(b_i_vec[b_j]) = t_dist(b_j);
        is_first_time[b_i_vec[b_j]] = false;
      }

    }
  }

  return tot;
}

// Compute nodal distance
double Nodal::GetNodalDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  if (path.GetMode()==NEW_MODE || is_first_time[b_i]) {
    switch(dist_type) {
      case 0:
        return NewtonRaphsonDistance(b_i,ref_b,other_b_i);
        break;
      case 1:
        return LineSearchDistance(b_i,ref_b,other_b_i);
        break;
      case 2:
        return MaxDistance(b_i,ref_b,other_b_i);
        break;
      case 3:
        return HybridDistance(b_i,ref_b,other_b_i);
        break;
      case 4:
        return HyperplaneDistance(b_i,ref_b,other_b_i);
        break;
    }
  } else
    return dist_c(b_i);

}

// First compute a single step of Newton-Raphson. Accept this distance if it is large (but smaller than maximum possible). Otherwise do a line search.
double Nodal::HybridDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  // Calculate determinant and gradiant
  SetRhoFGradRhoF(b_i,ref_b,other_b_i);
  if (rho_f(b_i) < 0.)
    return -1.;

  // Do single Newton-Raphson iteration
  double grad_mag(FieldVecMag(grad_rho_f,b_i));
  double dist_b_i = rho_f(b_i)/grad_mag;
  if (dist_b_i > 4.*sqrt(1./i_4_lambda_tau)) {
    double max_dist = MaxDistance(b_i,ref_b,other_b_i);
    return (dist_b_i < max_dist) ? dist_b_i : max_dist;
  } else
    return LineSearchDistance(b_i,ref_b,other_b_i);
}

// Return the slope and intercept of a line from 2 points
void Nodal::GetLine(vec<double> &p_a, vec<double> &p_b, double &m, double &b)
{
  if (p_a(1) > p_b(1))
    m = (p_a(1)-p_b(1))/(p_a(0)-p_b(0));
  else
    m = (p_b(1)-p_a(1))/(p_b(0)-p_a(0));
  b = p_a(1) - m*p_a(0);
}

// Return the distance to a hyperplane formed by the other particles
double Nodal::HyperplaneDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  double dist_b_i = 1./dist_tolerance;
  for (uint32_t p_i=0; p_i<n_part-1; ++p_i) {
    for (uint32_t p_j=p_i+1; p_j<n_part; ++p_j) {
      vec<double> r_tau_0 = path.GetR(path(species_i,p_i,b_i));
      vec<double> r_tau_1 = path.GetR(path(species_i,p_j,b_i));
      vec<double> r_star_0 = path.GetR(ref_b[p_i]);
      vec<double> r_star_1 = path.GetR(ref_b[p_j]);
      vec<double> r_tau_0_p(r_tau_0);
      r_tau_0_p(1) = r_tau_1(0) - (r_tau_0(1)-r_tau_1(1))*(r_star_0(1)-r_star_1(1))/(r_star_0(0)-r_star_1(0));
      double a, b(-1.), c;
      GetLine(r_tau_0_p, r_tau_1, a, c);
      dist_b_i = fabs(a*r_tau_0(0) + b*r_tau_0(1) + c)/sqrt(a*a + b*b);
      dist_b_i -= nearbyint(dist_b_i*path.iL)*path.L;
      dist_b_i = fabs(dist_b_i);
    }
  }
  return dist_b_i*M_SQRT2;
}

// Return the maximum distance to a node (coincidence point)
double Nodal::MaxDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  double dist_b_i = 1./dist_tolerance;
  for (uint32_t p_i=0; p_i<n_part-1; ++p_i)
    for (uint32_t p_j=p_i+1; p_j<n_part; ++p_j)
      dist_b_i = std::min(mag(path.Dr(path(species_i,p_i,b_i),path(species_i,p_j,b_i))), dist_b_i);
  dist_b_i *= M_SQRT1_2;
  return dist_b_i; // FIXME: Should I include ref distance?
}

// Perform a 1D line search in the direction of the gradient
double Nodal::LineSearchDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  // Calculate determinant and gradiant
  SetRhoFGradRhoF(b_i,ref_b,other_b_i);
  if (rho_f(b_i) < 0.)
    return -1.;

  // Set up temporary storage fields
  field<vec<double>> other_b_i_0(n_part);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    other_b_i_0(p_i) = path.GetR(other_b_i[p_i]);

  // Set initial gradient unit vector
  field<vec<double>> grad_rho_f_0(grad_rho_f.row(b_i));
  double grad_mag(FieldVecMag(grad_rho_f_0));
  double rho_f_0 = rho_f(b_i);
  double dist_0 = rho_f_0/grad_mag;
  double dist_0_i_grad_mag = dist_0/grad_mag;
  double abs_dist_0 = fabs(dist_0);

  // Find close point on other side of the node (along this distance of the gradient)
  bool done = false;
  double min_factor(0.), max_factor(0.5), max_dist(MaxDistance(b_i,ref_b,other_b_i));
  while ((rho_f_0*rho_f(b_i) > 0.) && (max_factor*dist_0 < max_dist)) {
    max_factor *= 2.;
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      path.GetR(other_b_i[p_i]) = other_b_i_0(p_i) - max_factor*dist_0_i_grad_mag*grad_rho_f_0(p_i);
    SetRhoF(b_i,ref_b,other_b_i);
  }

  // Perform a bisection search for the sign change
  double dist_b_i(max_dist);
  if (rho_f_0*rho_f(b_i) < 0.) {
    double try_factor(1.);
    while (((max_factor-min_factor)*abs_dist_0>dist_tolerance) && (min_factor*abs_dist_0<max_dist)) {
      try_factor = 0.5*(max_factor+min_factor);
      for (uint32_t p_i=0; p_i<n_part; ++p_i)
        path.GetR(other_b_i[p_i]) = other_b_i_0(p_i) - try_factor*dist_0_i_grad_mag*grad_rho_f_0(p_i);
      SetRhoF(b_i,ref_b,other_b_i);
      if (rho_f_0*rho_f(b_i) > 0.)
        min_factor = try_factor;
      else
        max_factor = try_factor;
    }
    if (min_factor*abs_dist_0 > max_dist)
      dist_b_i = max_dist;
    else
      dist_b_i = abs_dist_0*try_factor;
  }

  // Reset positions, determinant, and gradient
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    path.GetR(other_b_i[p_i]) = other_b_i_0(p_i);
  //grad_rho_f.row(b_i) = grad_rho_f_0;
  rho_f(b_i) = rho_f_0;

  return dist_b_i;
}

// Perform iterations of the Newton-Raphson root-finding algorithm to find the distance to the nearest node
double Nodal::NewtonRaphsonDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  // Calculate determinant and gradiant
  SetRhoFGradRhoF(b_i,ref_b,other_b_i);
  if (rho_f(b_i) < 0.)
    return -1.;

  // Set up temporary storage fields
  field<vec<double>> t_other_b_i(n_part), other_b_i_0(n_part);
  for (uint32_t p_i=0; p_i<n_part; ++p_i) {
    other_b_i_0(p_i) = path.GetR(other_b_i[p_i]);
    t_other_b_i(p_i) = path.GetR(other_b_i[p_i]);
  }

  // Set initial gradient unit vector
  field<vec<double>> grad_rho_f_0(grad_rho_f.row(b_i));
  double grad_mag(FieldVecMag(grad_rho_f_0));
  double rho_f_0 = rho_f(b_i);
  double dist_0 = rho_f_0/grad_mag;

  // Perform Newton-Raphson steps
  bool done = false;
  int i_dist_step = 0;
  while (!done && (i_dist_step < max_dist_steps)) {
    // Find zero crossing with current gradient
    double t_dist = 2.*rho_f(b_i)/grad_mag;
    do {
      t_dist *= 0.5;
      for (uint32_t p_i=0; p_i<n_part; ++p_i)
        path.GetR(other_b_i[p_i]) = t_other_b_i(p_i) - (t_dist/grad_mag)*grad_rho_f(b_i,p_i);
      SetRhoF(b_i,ref_b,other_b_i);
    } while (rho_f(b_i) < 0.);

    // Set new temporary path
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      t_other_b_i(p_i) = path.GetR(other_b_i[p_i]);

    // Check if tolerance is satisfied
    if (t_dist < dist_tolerance*dist_0)
      done = true;

    // Set new determinant and gradient
    SetRhoFGradRhoF(b_i,ref_b,other_b_i);
    grad_mag = FieldVecMag(grad_rho_f,b_i);

    // Iterate
    i_dist_step++;
  }

  // Compute distance
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    t_other_b_i(p_i) -= other_b_i_0(p_i);
  double dist_b_i(FieldVecMag(t_other_b_i));

  // Reset positions, determinant, and gradient
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    path.GetR(other_b_i[p_i]) = other_b_i_0(p_i);
  grad_rho_f.row(b_i) = grad_rho_f_0;
  rho_f(b_i) = rho_f_0;

  return dist_b_i;
}

// Form rho_f and grad_rho_f
void Nodal::SetRhoFGradRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
{
  // Set slice distance
  uint32_t abs_slice_diff = abs(path.bead_loop(b_i)-path.ref_bead);
  uint32_t min_slice_diff = std::min(abs_slice_diff, path.n_bead-abs_slice_diff);

  // Compute Gij matrix and its gradient
  mat<double> g_ij(n_part,n_part);
  field<vec<double>> dg_ij_dr(n_part,n_part);
  for (uint32_t p_i=0; p_i<n_part; ++p_i) {
    for (uint32_t p_j=0; p_j<n_part; ++p_j) {
      dg_ij_dr(p_i,p_j).zeros(path.n_d);
      g_ij(p_i,p_j) = GetGijDGijDr(path.Dr(ref_b[p_i], other_b_i[p_j]), min_slice_diff, dg_ij_dr(p_i,p_j));
    }
  }

  // Compute determinant with cofactors using Gauss-Jordan elimination
  rho_f(b_i) = det(g_ij);
  if (rho_f(b_i) < 0.) // Check sign
    return;

  // Compute gradient of determinant
  mat<double> i_g_ij(inv(g_ij));
  mat<double> cofactors_g_ij(rho_f(b_i)*trans(i_g_ij)); // TODO: Should just compute directly from adjugate
  for (uint32_t p_i=0; p_i<n_part; ++p_i) {
    grad_rho_f(b_i,p_i).zeros();
    for (uint32_t p_j=0; p_j<n_part; ++p_j)
      grad_rho_f(b_i,p_i) -= dg_ij_dr(p_j,p_i)*cofactors_g_ij(p_j,p_i); // FIXME: Should this be minue? Also, could be written in fewer lines.
  }

}

// Compute importance weight
double Nodal::ImportanceWeight()
{
  if (is_importance_weight) {
    std::vector<std::pair<uint32_t,uint32_t>> particles;
    particles.push_back(std::make_pair(species_i,0));
    path.SetMode(NEW_MODE);
    return GetAction(0, path.n_bead, particles, 0);
  } else
    return 1.;
}

// Accept new determinants and distances
void Nodal::Accept()
{
  if (use_nodal_distance) {
    for (uint32_t b_i=start_b; b_i<=end_b; ++b_i) {
      uint32_t t_b_i = path.bead_loop(b_i);
      rho_f_c(t_b_i) = rho_f(t_b_i);
      grad_rho_f_c.row(t_b_i) = grad_rho_f.row(t_b_i);
      dist_c(t_b_i) = dist(t_b_i);
    }
  }
}

// Reject new determinants and distances
void Nodal::Reject()
{
  if (use_nodal_distance) {
    for (uint32_t b_i=start_b; b_i<=end_b; ++b_i) {
      uint32_t t_b_i = path.bead_loop(b_i);
      rho_f(t_b_i) = rho_f_c(t_b_i);
      grad_rho_f.row(t_b_i) = grad_rho_f_c.row(t_b_i);
      dist(t_b_i) = dist_c(t_b_i);
    }
  }
}
