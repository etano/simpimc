#ifndef SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

/// Simple diagonal pair action
class BarePairAction : public PairAction
{
private:
  double k_cutoff; ///< Cutoff in k space
  double r_v_max; ///< Maximum distance possible in the pair potential spline
  double r_v_min; ///< Minimum distance possible in the pair potential spline
  double r_v_long_max; ///< Maximum distance possible in the long range pair potential spline
  double r_v_long_min; ///< Minimum distance possible in the long range pair potential spline
  double v_long_k_0; ///< Constant correction to the k components of the long range pair potential
  double v_long_r_0; ///< Constant correction to spatial long range pair potential
  NUBspline_1d_d *v_r_spline; ///< Spline of the pair potential
  NUBspline_1d_d *v_long_r_spline; ///< Spline of the long range pair potential
  vec<double> v_long_k; ///< Vector of k components of the long range pair potential
public:
  // Constructor
  BarePairAction(Path &path, Input &in, IO &out)
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
      v_long_k_0 *= 0.5*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      v_long_r_0 *= -0.5*species_a->GetNPart()*path.GetNBead();
    } else { // heterologous
      v_long_k_0 *= species_a->GetNPart()*species_b->GetNPart()*path.GetNBead();
      v_long_r_0 *= 0.;
    }

  }

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
        tot += CMag2(rho_k_a(b_i)(k_i),rho_k_b(b_i)(k_i))*v_long_k(k_i);

    if (species_a != species_b)
      tot *= 2.;

    return tot + v_long_k_0 + v_long_r_0;
  }

  /// Calculate the action
  virtual double CalcU(double r, double r_p, double s, const uint32_t level)
  {
    uint32_t skip = 1>>level;
    double level_tau = skip*path.GetTau();
    return level_tau*CalcV(r,r_p,level);
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
        tot += v_long_k(k_i)*CMag2(rho_k_a(species_a->bead_loop(b_i))(k_i),rho_k_b(species_b->bead_loop(b_i))(k_i));

    if (species_a != species_b)
      tot *= 2.;

    double level_tau = skip*path.GetTau();
    return level_tau*tot;
  }

  /// Calculate the beta derivative of the action
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level)
  {
    return CalcV(r,r_p,level);
  }

  /// Calculate the long ranged part of the beta derivative of the action
  virtual double CalcdUdBetaLong()
  {
    return CalcVLong();
  }

};

#endif // SIMPIMC_ACTIONS_BARE_PAIR_ACTION_CLASS_H_
