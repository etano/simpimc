#include "bisect_class.h"

// Create a spline for each possible slice_diff
void Bisect::SetupSpline()
{
  // Create splines
  uint32_t n_spline = floor(log2(path.n_bead));
  #pragma omp parallel for
  for (uint32_t spline_i=0; spline_i<n_spline; ++spline_i) {
    rho_free_splines.emplace_back(path.L, n_images, lambda, 0.5*path.tau*(1<<spline_i), false);
  }
}

void Bisect::Init(Input &in)
{
  n_level = in.GetAttribute<uint32_t>("n_level");
  uint32_t max_possible_level = floor(log2(path.n_bead));
  if (n_level > max_possible_level)
    std::cout << "Warning: n_level > max_possible_level!" << std::endl;
  if (path.pbc)
    n_images = in.GetAttribute<int>("n_images");
  else
    n_images = 0;

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
  Reset();
}

// Reset counters
void Bisect::Reset()
{
  if (!first_time && adaptive) {
    double accept_ratio = (double) n_accept / (double) n_attempt;
    if (accept_ratio < target_ratio && n_level > 1)
      n_level--;
    else if (1<<n_level < path.n_bead/2)
      n_level++;
    n_bisect_beads = 1<<n_level; // Number of beads in bisection
    i_4_lambda_tau_n_bisect_beads = i_4_lambda_tau/n_bisect_beads;
  }

  ref_accept = 0;
  ref_attempt = 0;

  Move::Reset();
}

// Accept current move
void Bisect::Accept()
{
  // Move Accepted, so copy new coordinates
  path.StoreR(affected_beads);
  path.StoreRhoKP(affected_beads);
  for (uint32_t b_i=bead0; b_i<bead1; ++b_i)
    path.StoreRhoK(b_i,species_i);

  // Call accept for each action
  for (auto& action: action_list)
    action->Accept();
}

// Reject current move
void Bisect::Reject()
{
  // Move rejected, so return old coordinates
  path.RestoreR(affected_beads);
  path.RestoreRhoKP(affected_beads);
  for (uint32_t b_i=bead0; b_i<bead1; ++b_i)
    path.RestoreRhoK(b_i,species_i);

  // Call reject for each action
  for (auto& action: action_list)
    action->Reject();
}

// Bisection Move
bool Bisect::Attempt()
{
  int p_i = rng.UnifRand(n_part) - 1;  // Pick particle at random
  bead0 = rng.UnifRand(path.n_bead) - 1;  // Pick first bead at random
  bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
  bool roll_over = bead1 > (path.n_bead-1);  // See if bisection overflows to next particle
  bool include_ref = path.species_list[species_i]->fixed_node &&
                    ((bead0<=path.ref_bead && bead1>=path.ref_bead) ||
                    (roll_over && path.bead_loop[bead1]>=path.ref_bead));
  if (include_ref)
    ref_attempt++;

  // Set up pointers
  std::shared_ptr<Bead> beadI(path(species_i,p_i,bead0));
  std::shared_ptr<Bead> beadF(beadI);
  affected_beads.clear();
  affected_beads.push_back(beadI);
  for (uint32_t i=0; i<n_bisect_beads; ++i) {
    beadF = beadF->next;
    affected_beads.push_back(beadF);
  }

  // Set which particles are affected by the move
  std::vector<std::pair<uint32_t,uint32_t>> particles;
  particles.push_back(std::make_pair(species_i,p_i));
  if (beadF->p != p_i)  // fixme: may be overkill
    particles.push_back(std::make_pair(species_i,beadF->p));

  // Perform the bisection (move exactly through kinetic action)
  std::shared_ptr<Bead> beadA, beadB, beadC;
  double prev_action_change = 0.;
  for (int level_i = n_level-1; level_i >= 0; level_i -= 1) {
    uint32_t skip = 1<<level_i;
    double level_tau = path.tau*skip;
    double sigma2 = lambda*level_tau;
    double sigma = sqrt(sigma2);

    double old_log_sample_prob = 0.;
    double new_log_sample_prob = 0.;
    beadA = beadI;
    while(beadA != beadF) {
      // Set beads
      beadB = beadA->NextB(skip);
      beadC = beadB->NextB(skip);

      // Old sampling
      path.SetMode(OLD_MODE);
      vec<double> delta_old(path.Dr(beadB, path.RBar(beadC, beadA)));
      old_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_old);

      // New sampling
      path.SetMode(NEW_MODE);
      vec<double> delta_new(path.n_d);
      rng.NormRand(delta_new, 0., sigma);
      path.PutInBox(delta_new);
      beadB->r = path.RBar(beadC, beadA) + delta_new;
      new_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_new);

      // Advance beads
      beadA = beadC;
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

    double log_sample_ratio = -new_log_sample_prob + old_log_sample_prob;
    double current_action_change = new_action - old_action;
    double log_accept_probablity = log_sample_ratio - current_action_change + prev_action_change;

    // Metropolis reject step
    if (log_accept_probablity < log(rng.UnifRand()))
      return 0;

    prev_action_change = current_action_change;
  }

  if (include_ref)
    ref_accept++;

  return 1;
}

void Bisect::Write()
{
  // Write
  if (first_time) {
    if (path.species_list[species_i]->fixed_node) {
      out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
      out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
    }
    if (adaptive)
      out.CreateExtendableDataSet("/Moves/"+name+"/", "n_level", n_level);
  } else {
    if (path.species_list[species_i]->fixed_node) {
      out.AppendDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
      out.AppendDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
    }
    if (adaptive)
      out.AppendDataSet("/Moves/"+name+"/", "n_level", n_level);
  }

  Move::Write();
}
