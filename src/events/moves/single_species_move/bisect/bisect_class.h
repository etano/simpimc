#ifndef SIMPIMC_MOVES_BISECT_CLASS_H_
#define SIMPIMC_MOVES_BISECT_CLASS_H_

#include "../single_species_move_class.h"
#include "../../../../actions/free_spline_class.h"

/// Performs a bisection move
class Bisect : public SingleSpeciesMove
{
protected:
  bool adaptive; ///< Whether or not to dynamically adjust the bisection level
  double i_4_lambda_tau_n_bisect_beads; ///< 1/(4\lambda\tau*n_bisect_beads) for affected species
  double target_ratio; ///< Target accepted ratio if using adaptive bisection level
  int n_images; ///< Number of periodic images used for the bisection sampling
  int ref_accept; ///< Counter tracking the number of bisections involving the reference slice are accepted
  int ref_attempt; ///< Counter tracking the number of bisections involving the reference slice are attempted
  uint32_t bead0; ///< First time slice in bisection
  uint32_t bead1; ///< Last time slice in bisection
  uint32_t n_bisect_beads; ///< Number of beads involved in a bisection of level n_level
  uint32_t n_level; ///< Number of levels in the bisection
  std::vector<FreeSpline> rho_free_splines; ///< Holds the splined action for every time slice
  std::vector<std::shared_ptr<Bead>> affected_beads; ///< Vector of the beads affected by the current bisection

  /// Accept the move
  virtual void Accept()
  {
    // Move Accepted, so copy new coordinates
    for (auto& b: affected_beads){
      b->StoreR();
      b->StoreRhoK();
    }
    for (uint32_t b_i=bead0; b_i<bead1; ++b_i)
      species->StoreRhoK(b_i);

    // Call accept for each action
    for (auto& action: action_list)
      action->Accept();
  }

  /// Attempt the move
  virtual bool Attempt()
  {
    int p_i = rng.UnifRand(species->GetNPart()) - 1;  // Pick particle at random
    bead0 = rng.UnifRand(species->GetNBead()) - 1;  // Pick first bead at random
    bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
    bool roll_over = bead1 > (species->GetNBead()-1);  // See if bisection overflows to next particle
    bool include_ref = species->IsFixedNode() &&
                      ((bead0<=species->GetRefBead() && bead1>=species->GetRefBead()) ||
                      (roll_over && species->bead_loop(bead1)>=species->GetRefBead()));
    if (include_ref)
      ref_attempt++;

    // Set up pointers
    std::shared_ptr<Bead> bead_i(species->GetBead(p_i,bead0));
    std::shared_ptr<Bead> bead_f(bead_i);
    affected_beads.clear();
    affected_beads.push_back(bead_i);
    for (uint32_t i=0; i<n_bisect_beads; ++i) {
      bead_f = bead_f->GetNextBead(1);
      affected_beads.push_back(bead_f);
    }

    // Set which particles are affected by the move
    std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> particles;
    particles.push_back(std::make_pair(species,p_i));
    if (bead_f->GetP() != p_i)  // TODO: may be overkill
      particles.push_back(std::make_pair(species,bead_f->GetP()));

    // Perform the bisection (move exactly through kinetic action)
    std::shared_ptr<Bead> bead_a, bead_b, bead_c;
    double prev_action_change = 0.;
    for (int level_i = n_level-1; level_i >= 0; level_i -= 1) {
      uint32_t skip = 1<<level_i;
      double level_tau = path.GetTau()*skip;
      double sigma2 = species->GetLambda()*level_tau;
      double sigma = sqrt(sigma2);

      double old_log_sample_prob = 0.;
      double new_log_sample_prob = 0.;
      bead_a = bead_i;
      while(bead_a != bead_f) {
        // Set beads
        bead_b = bead_a->GetNextBead(skip);
        bead_c = bead_b->GetNextBead(skip);

        // Old sampling
        path.SetMode(OLD_MODE);
        vec<double> delta_old(path.Dr(bead_b,path.RBar(bead_c,bead_a)));
        old_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_old);

        // New sampling
        path.SetMode(NEW_MODE);
        vec<double> delta_new(path.GetND());
        rng.NormRand(delta_new, 0., sigma);
        path.PutInBox(delta_new);
        bead_b->SetR(path.RBar(bead_c, bead_a) + delta_new);
        new_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_new);

        // Advance beads
        bead_a = bead_c;
      }

      // Calculate action change
      double old_action = 0.;
      double new_action = 0.;
      for (auto& action: action_list) {
        path.SetMode(OLD_MODE);
        old_action += action->GetAction(bead0,bead1,particles,level_i);
        path.SetMode(NEW_MODE);
        new_action += action->GetAction(bead0,bead1,particles,level_i);
      }

      // Metropolis reject step
      double log_sample_ratio = -new_log_sample_prob + old_log_sample_prob;
      double current_action_change = new_action - old_action;
      double log_accept_probability = log_sample_ratio - current_action_change + prev_action_change;
      if (log_accept_probability < log(rng.UnifRand()))
        return 0;

      prev_action_change = current_action_change;
    }

    if (include_ref)
      ref_accept++;

    return 1;
  }

  /// Rejects the move
  virtual void Reject()
  {
    // Move rejected, so return old coordinates
    for (auto& b: affected_beads){
      b->RestoreR();
      b->RestoreRhoK();
    }
    for (uint32_t b_i=bead0; b_i<bead1; ++b_i)
      species->RestoreRhoK(b_i);

    // Call reject for each action
    for (auto& action: action_list)
      action->Reject();
  }

  /// Resets the relevant counters
  virtual void Reset()
  {
    if (!first_time && adaptive) {
      double accept_ratio = (double) n_accept / (double) n_attempt;
      if (accept_ratio < target_ratio && n_level > 1)
        n_level--;
      else if (1<<n_level < species->GetNBead()/2)
        n_level++;
      n_bisect_beads = 1<<n_level; // Number of beads in bisection
      i_4_lambda_tau_n_bisect_beads = i_4_lambda_tau/n_bisect_beads;
    }

    ref_accept = 0;
    ref_attempt = 0;
  }

  /// Creates splined action for all time slices
  virtual void SetupSpline()
  {
    // Create splines
    uint32_t n_spline = floor(log2(species->GetNBead()));
    for (uint32_t spline_i=0; spline_i<n_spline; ++spline_i)
      rho_free_splines.emplace_back(path.GetL(), n_images, species->GetLambda(), 0.5*path.GetTau()*(1<<spline_i), false);
  }

public:
  /// Constructor reads in relevant inputs
  Bisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {
    n_level = in.GetAttribute<uint32_t>("n_level");
    uint32_t max_possible_level = floor(log2(species->GetNBead()));
    if (n_level > max_possible_level)
      std::cout << "Warning: n_level > max_possible_level!" << std::endl;
    n_images = in.GetAttribute<int>("n_images",0);

    // Adaptive bisection level
    adaptive = in.GetAttribute<bool>("adaptive",false);
    if (adaptive)
      target_ratio = in.GetAttribute<double>("target_ratio");

    // Compute constants
    n_bisect_beads = 1<<n_level; // Number of beads in bisection
    i_4_lambda_tau_n_bisect_beads = i_4_lambda_tau/n_bisect_beads;

    // Setup splines
    SetupSpline();

    // Reset counters
    Bisect::Reset();
    Move::Reset();
  }

  /// Writes relevant information about the move to the output file
  virtual void Write()
  {
    // Write
    if (first_time) {
      if (species->IsFixedNode()) {
        out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
        out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
      }
      if (adaptive)
        out.CreateExtendableDataSet("/Moves/"+name+"/", "n_level", n_level);
    } else {
      if (species->IsFixedNode()) {
        out.AppendDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
        out.AppendDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
      }
      if (adaptive)
        out.AppendDataSet("/Moves/"+name+"/", "n_level", n_level);
    }

    Bisect::Reset();
    Move::Write();
  }
};

#endif // SIMPIMC_MOVES_BISECT_CLASS_H_
