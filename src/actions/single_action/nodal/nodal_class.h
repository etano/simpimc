#ifndef SIMPIMC_ACTIONS_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_NODAL_CLASS_H_

#include "../single_action_class.h"

/// The nodal action class
class Nodal : public SingleAction
{
private:
  bool first_time_write; ///< Whether or not this is the first time writing to file
  uint32_t n_line_search; ///< Number of line searches performed
  uint32_t n_newton_raphson; ///< Number of complete newton raphson searches performed
protected:
  bool use_nodal_distance; ///< Whether or not to use a nodal distance
  double dist_tolerance; ///< Tolerance used when computing nodal distance
  int dist_type; ///< Index of type of distance measure used
  int max_dist_steps; ///< Maximum number of steps used to compute nodal distance
  uint32_t end_b; ///< Last bead affected by the action
  uint32_t start_b; ///< First bead affected by the action
  std::vector<bool> is_first_time; ///< Whether or not a specific time slice has been visisted
  vec<double> dist; ///< Vector of nodal distances for each time slice
  vec<double> dist_c; ///< Copy of vector of nodal distances for each time slice
  vec<double> rho_f; ///< Vector of nodal determinants for each time slice
  vec<double> rho_f_c; ///< Copy of vector of nodal determinants for each time slice
  field<vec<double>> grad_rho_f; ///< Field of nodal gradient vectors for each time slice
  field<vec<double>> grad_rho_f_c; ///< Copy of field of nodal gradient vectors for each time slice

  /// Compute the slope and y-intercept of a line given two points
  void GetLine(vec<double> &p_a, vec<double> &p_b, double &m, double &b)
  {
    if (p_a(1) > p_b(1))
      m = (p_a(1)-p_b(1))/(p_a(0)-p_b(0));
    else
      m = (p_b(1)-p_a(1))/(p_b(0)-p_a(0));
    b = p_a(1) - m*p_a(0);
  }

  /// Return the nodal distance for time slice b_i
  double GetNodalDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
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

  /// Compute the nodal distance by combining several methods
  double HybridDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
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

  /// Compute the nodal distance of time slice b_i by performing a 1D line search in the direction of the gradient
  double LineSearchDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    // Calculate determinant and gradiant
    SetRhoFGradRhoF(b_i,ref_b,other_b_i);
    if (rho_f(b_i) < 0.)
      return -1.;

    // Set up temporary storage fields
    field<vec<double>> other_b_i_0(species->GetNPart());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      other_b_i_0(p_i) = other_b_i[p_i]->GetR();

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
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b_i[p_i]->GetR() = other_b_i_0(p_i) - max_factor*dist_0_i_grad_mag*grad_rho_f_0(p_i);
      SetRhoF(b_i,ref_b,other_b_i);
    }

    // Perform a bisection search for the sign change
    double dist_b_i(max_dist);
    if (rho_f_0*rho_f(b_i) < 0.) {
      double try_factor(1.);
      while (((max_factor-min_factor)*abs_dist_0>dist_tolerance) && (min_factor*abs_dist_0<max_dist)) {
        try_factor = 0.5*(max_factor+min_factor);
        for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
          other_b_i[p_i]->GetR() = other_b_i_0(p_i) - try_factor*dist_0_i_grad_mag*grad_rho_f_0(p_i);
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

    // Reset positions and determinant
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      other_b_i[p_i]->GetR() = other_b_i_0(p_i);
    rho_f(b_i) = rho_f_0;

    // Track number of complete line search calls
    n_line_search++;

    return dist_b_i;
  }

  /// Compute the maximum nodal distance of time slice b_i
  double MaxDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    double dist_b_i = 1./dist_tolerance;
    for (uint32_t p_i=0; p_i<species->GetNPart()-1; ++p_i)
      for (uint32_t p_j=p_i+1; p_j<species->GetNPart(); ++p_j)
        dist_b_i = std::min(mag(path.Dr(species->GetBead(p_i,b_i),species->GetBead(p_j,b_i))), dist_b_i);
    dist_b_i *= M_SQRT1_2;
    return dist_b_i;
  }

  /// Compute the nodal distance of time slice b_i by finding the distance of a point to the hyperplane formed by the other particles
  double HyperplaneDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    double dist_b_i = 1./dist_tolerance;
    for (uint32_t p_i=0; p_i<species->GetNPart()-1; ++p_i) {
      for (uint32_t p_j=p_i+1; p_j<species->GetNPart(); ++p_j) {
        vec<double> r_tau_0 = species->GetBead(p_i,b_i)->GetR();
        vec<double> r_tau_1 = species->GetBead(p_j,b_i)->GetR();
        vec<double> r_star_0 = ref_b[p_i]->GetR();
        vec<double> r_star_1 = ref_b[p_j]->GetR();
        vec<double> r_tau_0_p(r_tau_0);
        r_tau_0_p(1) = r_tau_1(0) - (r_tau_0(1)-r_tau_1(1))*(r_star_0(1)-r_star_1(1))/(r_star_0(0)-r_star_1(0));
        double a, b(-1.), c;
        GetLine(r_tau_0_p, r_tau_1, a, c);
        dist_b_i = fabs(a*r_tau_0(0) + b*r_tau_0(1) + c)/sqrt(a*a + b*b);
        dist_b_i -= nearbyint(dist_b_i*path.GetInverseL())*path.GetL();
        dist_b_i = fabs(dist_b_i);
      }
    }
    return dist_b_i*M_SQRT2;
  }

  /// Compute the nodal distance of time slice b_i by performing a newton raphson search
  double NewtonRaphsonDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    // Calculate determinant and gradiant
    SetRhoFGradRhoF(b_i,ref_b,other_b_i);
    if (rho_f(b_i) < 0.)
      return -1.;
    if (std::isnan(grad_rho_f(b_i,0)(0)))
      return LineSearchDistance(b_i,ref_b,other_b_i);

    // Set up temporary storage fields
    field<vec<double>> t_other_b_i(species->GetNPart()), other_b_i_0(species->GetNPart());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
      other_b_i_0(p_i) = other_b_i[p_i]->GetR();
      t_other_b_i(p_i) = other_b_i[p_i]->GetR();
    }

    // Set initial gradient unit vector
    field<vec<double>> grad_rho_f_0(grad_rho_f.row(b_i));
    double grad_mag(FieldVecMag(grad_rho_f_0));
    if (grad_mag < dist_tolerance)
      return LineSearchDistance(b_i,ref_b,other_b_i);
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
        for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
          other_b_i[p_i]->GetR() = t_other_b_i(p_i) - (t_dist/grad_mag)*grad_rho_f(b_i,p_i);
        }
        SetRhoF(b_i,ref_b,other_b_i);
      } while (rho_f(b_i) < 0.);

      // Set new temporary path
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        t_other_b_i(p_i) = other_b_i[p_i]->GetR();

      // Check if tolerance is satisfied
      if (t_dist < dist_tolerance*dist_0)
        done = true;

      // Set new determinant and gradient
      SetRhoFGradRhoF(b_i,ref_b,other_b_i);
      if (!std::isnan(grad_rho_f(b_i,0)(0))) {
        grad_mag = FieldVecMag(grad_rho_f,b_i);
        if (grad_mag > 0.) {
          // Iterate and continue
          i_dist_step++;
          continue;
        }
      }

      // Reset positions, determinant, and gradient
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b_i[p_i]->GetR() = other_b_i_0(p_i);
      grad_rho_f.row(b_i) = grad_rho_f_0;
      rho_f(b_i) = rho_f_0;
      return LineSearchDistance(b_i,ref_b,other_b_i);

    }

    // Compute distance
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
      t_other_b_i(p_i) -= other_b_i_0(p_i);
      path.PutInBox(t_other_b_i(p_i));
    }
    double dist_b_i(FieldVecMag(t_other_b_i));

    // Reset positions, determinant, and gradient
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      other_b_i[p_i]->GetR() = other_b_i_0(p_i);
    grad_rho_f.row(b_i) = grad_rho_f_0;
    rho_f(b_i) = rho_f_0;

    // Track number of complete Newton Raphson calls
    n_newton_raphson++;

    return dist_b_i;
  }

  /// Compute the determinant nodal action for time slice b_i
  void SetRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    // Set slice distance
    uint32_t abs_slice_diff = abs(species->bead_loop(b_i)-species->GetRefBead());
    uint32_t min_slice_diff = std::min(abs_slice_diff, species->GetNBead()-abs_slice_diff);

    // Compute Gij matrix
    mat<double> g_ij(species->GetNPart(),species->GetNPart());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
      for (uint32_t p_j=0; p_j<species->GetNPart(); ++p_j) {
        g_ij(p_i,p_j) = GetGij(ref_b[p_i], other_b_i[p_j], min_slice_diff); // TODO: Fast updates
      }
    }

    // Compute determinant
    rho_f(b_i) = det(g_ij);

    return;
  }

  /// Compute the gradient of the nodal action for time slice b_i
  void SetRhoFGradRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i)
  {
    // Set slice distance
    uint32_t abs_slice_diff = abs(species->bead_loop(b_i)-species->GetRefBead());
    uint32_t min_slice_diff = std::min(abs_slice_diff, species->GetNBead()-abs_slice_diff);

    // Compute Gij matrix and its gradient
    mat<double> g_ij(species->GetNPart(),species->GetNPart());
    field<vec<double>> dg_ij_dr(species->GetNPart(),species->GetNPart());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
      for (uint32_t p_j=0; p_j<species->GetNPart(); ++p_j) {
        dg_ij_dr(p_i,p_j).zeros(path.GetND());
        g_ij(p_i,p_j) = GetGijDGijDr(ref_b[p_i], other_b_i[p_j], min_slice_diff, dg_ij_dr(p_i,p_j));
      }
    }

    // Compute determinant with cofactors using Gauss-Jordan elimination
    rho_f(b_i) = det(g_ij);
    if (rho_f(b_i) < 0.) // Check sign
      return;

    // Create cofactor matrix
    mat<double> i_g_ij;
    try {
      i_g_ij = inv(g_ij);
    } catch (const std::runtime_error& error) {
      rho_f(b_i) = -1.;
      grad_rho_f(b_i,0) = sqrt(-1.);
      return;
    }

    mat<double> cofactors_g_ij(rho_f(b_i)*trans(i_g_ij)); // TODO: Should just compute directly from adjugate

    // Compute gradient of determinant
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
      grad_rho_f(b_i,p_i).zeros();
      for (uint32_t p_j=0; p_j<species->GetNPart(); ++p_j)
        grad_rho_f(b_i,p_i) -= dg_ij_dr(p_j,p_i)*cofactors_g_ij(p_j,p_i);
    }

  }

  /// Returns the value of g_ij
  virtual double GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff) = 0;

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double> &dgij_dr) = 0;

  /// Compute the nodal action without using the nodal distance
  double SimpleAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all)
  {
    // Compute action
    double tot = 0.;
      for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
        if (b_i_vec[b_i] != species->GetRefBead())  {
          SetRhoF(b_i_vec[b_i],ref_b,other_b[b_i]);
          if (rho_f(b_i_vec[b_i]) < 0.) {
            tot += 1.e100;
            break;
          }
        }
      }

    return tot;
  }

  /// Compute the nodal action by using a nodal distance measure
  double DistanceAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all)
  {
    // Set nodal distances
    double tot = 0.;
    vec<double> t_dist(n_bead_in_move);
      for (uint32_t b_i=0; b_i<n_bead_in_move; ++b_i) {
        if (b_i_vec[b_i] != species->GetRefBead())  {
          t_dist(b_i) = GetNodalDistance(b_i_vec[b_i],ref_b,other_b[b_i]);
          if (t_dist(b_i) < 0.) {
            tot += 1.e100;
            break;
          }
        }
      }

    // Compute action from nodal distance
    if (tot == 0.) {
      uint32_t skip = b_i_vec[1] - b_i_vec[0];
      double i_lambda_level_tau = 1./(species->GetLambda()*skip*path.GetTau());
      for (uint32_t b_i=0; b_i<n_bead_in_move-1; ++b_i) {

        // Set slice variables
        uint32_t b_j = b_i + 1;
        bool b_i_is_ref = b_i_vec[b_i] == species->GetRefBead();
        bool b_j_is_ref = b_i_vec[b_j] == species->GetRefBead();

        // Compute action
        if (!b_i_is_ref && (t_dist(b_i)<0.)) {
          tot = 1.e100;
          break;
        } else if (!b_j_is_ref && (t_dist(b_j)<0.)) {
          tot = 1.e100;
          break;
        } else if (b_i_is_ref || (t_dist(b_i)==0.))
          tot -= log1p(-exp(-t_dist(b_j)*t_dist(b_j)*i_lambda_level_tau));
        else if (b_j_is_ref || (t_dist(b_j)==0.))
          tot -= log1p(-exp(-t_dist(b_i)*t_dist(b_i)*i_lambda_level_tau));
        else
          tot -= log1p(-exp(-t_dist(b_i)*t_dist(b_j)*i_lambda_level_tau));

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

  /// Creates splined action for all time slices and parameter sets
  virtual void SetupSpline() = 0;

  /// Compute magnitude from field of vectors
  template<typename T>
  T FieldVecMag(const field<vec<T>> &field_vec)
  {
    T tot = 0.;
    for (uint32_t i=0; i<field_vec.size(); ++i)
      tot += dot(field_vec(i),field_vec(i));
    return sqrt(tot);
  }

  /// Compute magnitude from field of vectors
  template<typename T>
  T FieldVecMag(const field<vec<T>> &field_vec, const uint32_t row)
  {
    T tot = 0.;
    for (uint32_t i=0; i<field_vec.n_cols; ++i)
      tot += dot(field_vec(row,i),field_vec(row,i));
    return sqrt(tot);
  }

  /// Test if nodes are broken
  bool TestNodes()
  {
    // Test initial configuration
    is_first_time.resize(species->GetNBead());
    for (uint32_t b_i=0; b_i<species->GetNBead(); ++b_i)
      is_first_time[b_i] = true;
    std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> particles;
    particles.push_back(std::make_pair(species,0));
    bool init_good = true;
    ModeType t_mode(path.GetMode());
    path.SetMode(NEW_MODE);
    if (GetAction(0, species->GetNBead(), particles, 0) >= 1.e100) {
      std::cout << "Warning: initializing with broken nodes!" << std::endl;
      init_good = false;
    }
    start_b = 0;
    end_b = species->GetNBead();
    Accept();
    path.SetMode(t_mode);

    // Reset things
    first_time_write = true;
    n_newton_raphson = 0;
    n_line_search = 0;

    return init_good;
  }

public:
  /// Constructor only calls parent Action class constructor
  Nodal(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {
    // Read in things
    is_importance_weight = in.GetAttribute<bool>("is_importance_weight",false);

    // Nodal distance things
    use_nodal_distance = in.GetAttribute<bool>("use_nodal_distance",false);
    out.Write("Actions/"+name+"/use_nodal_distance", use_nodal_distance);
    if (use_nodal_distance) {
      std::string dist_type_name = in.GetAttribute<std::string>("dist_type");
      std::cout << "Setting up " << dist_type_name << " distance measure for nodal action..." << std::endl;
      if (dist_type_name == "NewtonRaphson") {
        dist_type = 0;
        max_dist_steps = in.GetAttribute<int>("max_dist_steps",10);
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
      out.Write("Actions/"+name+"/dist_type", dist_type);
      dist_tolerance = in.GetAttribute<double>("tolerance",1.e-4);
      out.Write("Actions/"+name+"/dist_tolerance", dist_tolerance);
      dist.zeros(species->GetNBead());
      dist_c.zeros(species->GetNBead());
    }

    // Set up determinants
    rho_f.zeros(species->GetNBead());
    rho_f_c.zeros(species->GetNBead());
    grad_rho_f.set_size(species->GetNBead(),species->GetNPart());
    grad_rho_f_c.set_size(species->GetNBead(),species->GetNPart());
    for (uint32_t b_i=0; b_i<species->GetNBead(); ++b_i) {
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
        grad_rho_f(b_i,p_i).zeros(path.GetND());
        grad_rho_f_c(b_i,p_i).zeros(path.GetND());
      }
    }

  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta()
  {
    // Zero if not using the nodal distance or using it as an importance weight
    if (!use_nodal_distance)
      return 0.;

    // Initialize other beads
    std::vector<std::shared_ptr<Bead>> ref_b(species->GetNPart());
    std::vector<std::vector<std::shared_ptr<Bead>>> other_b(species->GetNBead());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      ref_b[p_i] = species->GetBead(p_i,species->GetRefBead());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      other_b[0].push_back(ref_b[p_i]->GetPrevBead(species->GetRefBead()));

    std::vector<uint32_t> b_i_vec;
    b_i_vec.push_back(other_b[0][0]->GetB());
    for (uint32_t b_i=1; b_i<species->GetNBead(); ++b_i) {
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b[b_i].push_back(other_b[b_i-1][p_i]->GetNextBead(1));
      b_i_vec.push_back(other_b[b_i][0]->GetB());
    }

    // Compute action from nodal distance
    double tot = 0.;
    double i_lambda_tau = 1./(species->GetLambda()*path.GetTau());
    for (uint32_t b_i=0; b_i<species->GetNBead(); ++b_i) {
      // Set slice variables
      uint32_t b_j = species->bead_loop(b_i + 1);
      bool b_i_is_ref(b_i == species->GetRefBead());
      bool b_j_is_ref(b_j == species->GetRefBead());

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
        tot += dist2_i_lambda_tau/(path.GetTau()*exp_dist2_i_lambda_tau_m1);
    }

    return tot/species->GetNBead();
  }

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level)
  {
    // No need to check old action if not using the nodal distance
    if (path.GetMode()==OLD_MODE && !use_nodal_distance)
      return 0.;

    // Decide whether or not to check the node
    if (level > max_level || !species->IsFermi())
      return 0.;
    bool check_node = false;
    for (auto& p: particles) {
      if (p.first == species) {
        check_node = true;
        break;
      }
    }
    if (!check_node)
      return 0.;

    // See if ref slice included
    bool check_all = false;
    if (b1 < species->GetNBead())
      check_all = ((b0 <= species->GetRefBead()) && (b1 >= species->GetRefBead()));
    else
      check_all = (species->bead_loop(b1) >= species->GetRefBead());

    // Set start and end beads
    if (check_all) {
      start_b = 0;
      end_b = species->GetNBead()-1;
    } else {
      start_b = b0;
      end_b = b1;
    }

    // Set up bead index vector
    uint32_t skip = 1<<level;
    int n_bead_in_move = 1 + (end_b - start_b)/skip;
    std::vector<uint32_t> b_i_vec;
    for (uint32_t b_i=start_b; b_i<=end_b; b_i+=skip)
      b_i_vec.push_back(species->bead_loop(b_i));

    // Initialize other beads
    std::vector<std::shared_ptr<Bead>> ref_b(species->GetNPart());
    std::vector<std::vector<std::shared_ptr<Bead>>> other_b(n_bead_in_move);
    int slice_diff_0 = species->bead_loop(start_b) - species->GetRefBead();
    int abs_slice_diff_0 = abs(slice_diff_0);
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      ref_b[p_i] = species->GetBead(p_i,species->GetRefBead());
    if (slice_diff_0 >= 0) {
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b[0].push_back(ref_b[p_i]->GetNextBead(abs_slice_diff_0));
    } else {
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b[0].push_back(ref_b[p_i]->GetPrevBead(abs_slice_diff_0));
    }
    for (uint32_t b_i=1; b_i<n_bead_in_move; ++b_i)
      for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
        other_b[b_i].push_back(other_b[b_i-1][p_i]->GetNextBead(skip));

    // Decide which type of action to compute
    if (use_nodal_distance)
      return DistanceAction(b_i_vec, ref_b, other_b, n_bead_in_move, check_all);
    else
      return SimpleAction(b_i_vec, ref_b, other_b, n_bead_in_move, check_all);
  }

  /// Returns the importance weight of the action for the whole path
  virtual double ImportanceWeight()
  {
    if (is_importance_weight) {
      std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> particles;
      particles.push_back(std::make_pair(species,0));
      path.SetMode(NEW_MODE);
      return GetAction(0, species->GetNBead(), particles, 0);
    } else
      return 1.;
  }

  /// Writes information about the action
  virtual void Write()
  {
    // Write out and reset counters if performing newton raphson
    if (use_nodal_distance && dist_type==0) {
      double fraction_newton_raphson = double(n_newton_raphson)/double(n_line_search+n_newton_raphson);
      std::string prefix = "Actions/"+name+"/";
      if (first_time_write) {
        out.CreateGroup(prefix+"fraction_newton_raphson");
        out.CreateExtendableDataSet("/"+prefix+"fraction_newton_raphson/", "x", fraction_newton_raphson);
        std::string data_type = "scalar";
        out.Write(prefix+"fraction_newton_raphson/data_type",data_type);
        first_time_write = false;
      } else {
        out.AppendDataSet("/"+prefix+"fraction_newton_raphson/", "x", fraction_newton_raphson);
      }
      n_newton_raphson = 0;
      n_line_search = 0;
    }
  }

  /// Accepts relevant information about the action
  virtual void Accept()
  {
    if (use_nodal_distance) {
      for (uint32_t b_i=start_b; b_i<=end_b; ++b_i) {
        uint32_t t_b_i = species->bead_loop(b_i);
        rho_f_c(t_b_i) = rho_f(t_b_i);
        grad_rho_f_c.row(t_b_i) = grad_rho_f.row(t_b_i);
        dist_c(t_b_i) = dist(t_b_i);
      }
    }
  }

  /// Rejects relevant information about the action
  virtual void Reject()
  {
    if (use_nodal_distance) {
      for (uint32_t b_i=start_b; b_i<=end_b; ++b_i) {
        uint32_t t_b_i = species->bead_loop(b_i);
        rho_f(t_b_i) = rho_f_c(t_b_i);
        grad_rho_f.row(t_b_i) = grad_rho_f_c.row(t_b_i);
        dist(t_b_i) = dist_c(t_b_i);
      }
    }
  }

};

#endif // SIMPIMC_ACTIONS_NODAL_CLASS_H_
