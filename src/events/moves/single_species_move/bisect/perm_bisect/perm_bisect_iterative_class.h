#ifndef SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_
#define SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_

#include "perm_bisect_class.h"

/// Permuting bisection move that builds the permutation iteratively
class PermBisectIterative : public PermBisect
{
private:
  /// Update the permutation table
  void UpdatePermTable()
  {
    // Set initial and final beads
    field< std::shared_ptr<Bead>> b0(n_part), b1(n_part);
    for (uint32_t p_i=0; p_i<n_part; p_i++) {
      b0(p_i) = path(species_i,p_i,bead0);
      b1(p_i) = b0(p_i)->NextB(n_bisect_beads);
    }

    // Construct t table
    for (uint32_t i=0; i<n_part; i++) {
      for (uint32_t j=0; j<n_part; j++) {
        vec<double> dr_ij(path.Dr(b0(i), b1(j)));
        double exponent = (-dot(dr_ij,dr_ij))*i_4_lambda_tau_n_bisect_beads;
        if (exponent > log_epsilon)
          t(i,j) = exp(exponent);
        else
          t(i,j) = 0.;
      }
    }

  }

  /// Select which cycle to use in the permutation
  bool SelectCycleIterative(Cycle &c)
  {
    // Update t
    UpdatePermTable();
    mat<double> t_c = t;

    // Choose particles
    uint32_t p0 = rng.UnifRand(n_part) - 1;  // Pick first particle at random
    uint32_t p = p0;
    std::vector<uint32_t> ps;
    uint32_t n_perm = 0;
    do {
      // Add particle to ps
      ps.push_back(p);
      n_perm = ps.size();

      // Make sure returning to previous particles is not an option
      for (uint32_t i=0; i<n_perm; ++i)
        t_c(p,ps[i]) = 0.;

      // Allow cycle to close
      if (n_perm > 0)
        t_c(p,p0) = t(p,p0);

      // Disallow odd permutations (even number of particles) for fixed-node fermions
      bool is_fixed_node_odd_perm = path.species_list[species_i]->fermi && path.species_list[species_i]->fixed_node && !(n_perm%2);
      if (is_fixed_node_odd_perm)
        t_c(p,p0) = 0.;

      // Calculate row total
      double Q_p = 0.;
      double Q_p_c = 0.;
      for (uint32_t i=0; i<n_part; ++i) {
        Q_p += t(p,i);
        Q_p_c += t_c(p,i);
      }

      // Decide whether or not to continue
      if (Q_p_c/Q_p < rng.UnifRand()) {
        c.type = n_perm - 1;
        return 0;
      }

      // Select next particle with bisective search
      double x = rng.UnifRand();
      double t_Q = 0.;
      for (uint32_t i=0; i<n_part; ++i) { // TODO: not doing bisection
        t_Q += t_c(p,i)/Q_p_c;
        if (t_Q > x) {
          p = i;
          break;
        }
      }

    } while (p != p0);

    // Set weight
    c.weight = 1.;
    for (uint32_t i=0; i<n_perm-1; ++i)
      c.weight *= t(ps[i],ps[i+1])/t(ps[i],ps[i]);
    c.weight *= t(ps[n_perm-1],ps[0])/t(ps[n_perm-1],ps[n_perm-1]);

    // Set particles
    c.part.set_size(n_perm);
    for (uint32_t i=0; i<n_perm; ++i)
      c.part(i) = ps[i];
    c.type = n_perm-1;

    // Set perms
    c.perm.set_size(n_perm);
    for (uint32_t i=0; i<n_perm-1; ++i)
      c.perm(i) = i+1;
    c.perm(n_perm-1) = 0;
    c.i_perm.set_size(n_perm);
    c.i_perm(0) = n_perm-1;
    for (uint32_t i=1; i<n_perm; ++i)
      c.i_perm(i) = i-1;

    return 1;
  }

  /// Attempt the move
  virtual bool Attempt()
  {
    bead0 = rng.UnifRand(path.n_bead) - 1;  // Pick first bead at random
    bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
    bool roll_over = bead1 > (path.n_bead-1);  // See if bisection overflows to next particle
    bool include_ref = path.species_list[species_i]->fixed_node &&
                      ((bead0<=path.ref_bead && bead1>=path.ref_bead) ||
                      (roll_over && path.bead_loop[bead1]>=path.ref_bead));
    if (include_ref)
      ref_attempt++;

    // Set up permutation
    Cycle c;
    n_perm_part = 0; // reset to indicate if bisection is attempted or not
    if (!SelectCycleIterative(c)) {
      n_perm_part = 0;
      return 0; // do not attempt bisection since permutation not accepted
    }
    n_perm_part = c.type+1;

    // Set up pointers
    std::vector<std::pair<uint32_t,uint32_t>> particles;
    field< std::shared_ptr<Bead>> bead_i(n_perm_part), bead_fm1(n_perm_part), bead_f(n_perm_part);
    for (uint32_t i=0; i<n_perm_part; i++) {
      bead_i(i) = path(species_i,c.part(i),bead0);
      bead_fm1(i) = bead_i(i)->NextB(n_bisect_beads-1);
      bead_f(i) = bead_fm1(i)->next;
      particles.push_back(std::make_pair(species_i,c.part(i)));
    }
    if (roll_over) {
      for (uint32_t i=0; i<n_perm_part; ++i)
        particles.push_back(std::make_pair(species_i,bead_f(i)->p));
      sort(particles.begin(), particles.end());
      particles.erase(unique(particles.begin(), particles.end()), particles.end());
    }

    // Permute particles
    PermuteBeads(bead_fm1, bead_f, c);

    // Note affected beads
    field< std::shared_ptr<Bead>> bead_a(n_perm_part);
    affected_beads.clear();
    for (uint32_t i=0; i<n_perm_part; i++) {
      for(bead_a(i) = bead_i(i)->next; bead_a(i) != bead_f(i); bead_a(i) = bead_a(i)->next)
        affected_beads.push_back(bead_a(i));
    }

    // Perform the bisection (move exactly through kinetic action)
    field< std::shared_ptr<Bead>> bead_b(n_perm_part), bead_c(n_perm_part);
    double prev_action_change = -log(c.weight);
    for (int level_i = n_level-1; level_i >= 0; level_i -= 1) {

      // Level specific quantities
      uint32_t skip = 1<<level_i;
      double sigma = sqrt(path.tau*skip*lambda);

      // Calculate sampling probability
      double old_log_sample_prob = 0.;
      double new_log_sample_prob = 0.;
      for (uint32_t i=0; i<n_perm_part; i++) {
        bead_a(i) = bead_i(i);
        while(bead_a(i) != bead_f(i)) {
          // Old sampling
          path.SetMode(OLD_MODE);
          bead_b(i) = path.GetNextBead(bead_a(i),skip);
          bead_c(i) = path.GetNextBead(bead_b(i),skip);
          vec<double> delta_old(0.5*(path.Dr(bead_b(i),bead_c(i)) + path.Dr(bead_b(i),bead_a(i))));
          path.PutInBox(delta_old);
          old_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_old);

          // New sampling
          path.SetMode(NEW_MODE);
          bead_b(i) = path.GetNextBead(bead_a(i),skip);
          bead_c(i) = path.GetNextBead(bead_b(i),skip);
          vec<double> delta_new(path.n_d);
          rng.NormRand(delta_new, 0., sigma);
          path.PutInBox(delta_new);
          bead_b(i)->r = path.RBar(bead_c(i),bead_a(i)) + delta_new;
          new_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_new);

          bead_a(i) = bead_c(i);
        }
      }

      // Calculate action change
      double old_action = 0.;
      double new_action = 0.;
      for (auto& action: action_list) {
        // Old action
        path.SetMode(OLD_MODE);
        old_action += action->GetAction(bead0, bead1, particles, level_i);

        // New action
        path.SetMode(NEW_MODE);
        new_action += action->GetAction(bead0, bead1, particles, level_i);
      }

      // Calculate acceptance ratio
      double log_sample_ratio = -new_log_sample_prob + old_log_sample_prob;
      double current_action_change = new_action - old_action;
      double log_accept_probability = log_sample_ratio - current_action_change + prev_action_change;

      // Metropolis step
      if (log_accept_probability < log(rng.UnifRand()))
        return 0;

      prev_action_change = current_action_change;
    }

    if (include_ref)
      ref_accept++;

    return 1;
  }
public:
  /// Constructor instantiates parent class and calls Init
  PermBisectIterative(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : PermBisect(path, rng, action_list, in, out)
  {}

};

#endif // SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_
