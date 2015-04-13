#ifndef SIMPIMC_ACTIONS_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_NODAL_CLASS_H_

#include "single_action_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

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
  void GetLine(vec<double> &p_a, vec<double> &p_b, double &m, double &b);

  /// Return the nodal distance for time slice b_i
  double GetNodalDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the nodal distance by combining several methods
  double HybridDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the nodal distance of time slice b_i by performing a 1D line search in the direction of the gradient
  double LineSearchDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the maximum nodal distance of time slice b_i
  double MaxDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the nodal distance of time slice b_i by finding the distance of a point to the hyperplane formed by the other particles
  double HyperplaneDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the nodal distance of time slice b_i by performing a newton raphson search
  double NewtonRaphsonDistance(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the determinant nodal action for time slice b_i
  void SetRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Compute the gradient of the nodal action for time slice b_i
  void SetRhoFGradRhoF(const int b_i, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::shared_ptr<Bead>> &other_b_i);

  /// Returns the value of g_ij
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff) = 0;

  /// Returns the spatial derivative of g_ij
  virtual double GetGijDGijDr(const vec<double> &r, const uint32_t slice_diff, vec<double> &dgij_dr) = 0;

  /// Compute the nodal action without using the nodal distance
  double SimpleAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all);

  /// Compute the nodal action by using a nodal distance measure
  double DistanceAction(const std::vector<uint32_t> &b_i_vec, const std::vector<std::shared_ptr<Bead>> &ref_b, const std::vector<std::vector<std::shared_ptr<Bead>>> &other_b, const int n_bead_in_move, const bool check_all);

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

public:
  /// Constructor only calls parent Action class constructor
  Nodal(Path &path, Input &in, IO &out)
    : SingleAction(path,in,out)
  {}

  /// Initialize the action
  virtual void Init(Input &in);

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta();

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Returns the importance weight of the action for the whole path
  virtual double ImportanceWeight();

  /// Writes information about the action
  virtual void Write();

  /// Accepts relevant information about the action
  virtual void Accept();

  /// Rejects relevant information about the action
  virtual void Reject();
};

#endif // SIMPIMC_ACTIONS_NODAL_CLASS_H_
