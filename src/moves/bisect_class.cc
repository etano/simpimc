#include "bisect_class.h"

void Bisect::Init(Input &in)
{
  n_level = in.GetAttribute<uint>("n_level");
  uint max_possible_level = floor(log2(path.n_bead));
  if (n_level > max_possible_level)
    std::cout << "Warning: n_level > max_possible_level!" << std::endl;
  if (path.pbc)
    n_images = in.GetAttribute<int>("n_images");
  else
    n_images = 0;
  species = in.GetAttribute<std::string>("species");
  path.GetSpeciesInfo(species,species_i);

  // Generate action list
  std::vector<std::string> species_list;
  species_list.push_back(species);
  GenerateActionList(species_list);

  // Adaptive bisection level
  adaptive = in.GetAttribute<bool>("adaptive",false);
  if (adaptive)
    target_ratio = in.GetAttribute<double>("target_ratio");

  // Compute constants
  n_bisect_beads = 1<<n_level; // Number of beads in bisection
  lambda = path.species_list[species_i]->lambda;
  i4_lambda_tau_n_bisect_beads = 1./(4.*lambda*path.tau*n_bisect_beads);
}

// Reset counters
void Bisect::Reset()
{
  if (adaptive) {
    double accept_ratio = (double) n_accept / (double) n_attempt;
    if (accept_ratio < target_ratio && n_level > 1)
      n_level--;
    else if (1<<n_level < path.n_bead/2)
      n_level++;
    n_bisect_beads = 1<<n_level; // Number of beads in bisection
    lambda = path.species_list[species_i]->lambda;
    i4_lambda_tau_n_bisect_beads = 1./(4.*lambda*path.tau*n_bisect_beads);
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
  for (uint b_i=bead0+1; b_i<bead1; ++b_i)
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
  for (uint b_i=bead0+1; b_i<bead1; ++b_i)
    path.RestoreRhoK(b_i,species_i);

  // Call reject for each action
  for (auto& action: action_list)
    action->Reject();
}

// Bisection Move
bool Bisect::Attempt()
{
  uint p_i = rng.UnifRand(path.species_list[species_i]->n_part) - 1;  // Pick particle at random
  bead0 = rng.UnifRand(path.n_bead) - 1;  // Pick first bead at random
  bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
  roll_over = bead1 > (path.n_bead-1);  // See if bisection overflows to next particle
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
  for (uint i=0; i<n_bisect_beads; ++i) {
    beadF = beadF->next;
    affected_beads.push_back(beadF);
  }

  // Set which particles are affected by the move
  std::vector<std::pair<uint,uint>> particles;
  particles.push_back(std::make_pair(species_i,p_i));
  if (beadF->p != p_i)  // fixme: may be overkill
    particles.push_back(std::make_pair(species_i,beadF->p));

  // Perform the bisection (move exactly through kinetic action)
  std::shared_ptr<Bead> beadA, beadB, beadC;
  double prev_action_change = 0.;
  double prefactor_sample_prob = 0.;
  for (int level_i = n_level-1; level_i >= 0; level_i -= 1) {
    uint skip = 1<<level_i;
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
      path.SetMode(0);
      vec<double> r_bar_old(path.RBar(beadC, beadA));
      vec<double> delta_old(path.Dr(beadB, r_bar_old));

      // New sampling
      path.SetMode(1);
      vec<double> r_bar_new(path.RBar(beadC, beadA));
      vec<double> delta_new(path.n_d);
      rng.NormRand(delta_new, 0., sigma);
      path.PutInBox(delta_new);
      beadB->r = r_bar_new + delta_new;

      // Get sampling probs
      double gauss_prod_old = 1.;
      double gauss_prod_new = 1.;
      for (uint iD=0; iD<path.n_d; iD++) {
        double gauss_sum_old = 0.;
        double gauss_sum_new = 0.;
        for (int image=-n_images; image<=n_images; image++) {
          double dist_old = delta_old(iD) + (double)image*path.L;
          double dist_new = delta_new(iD) + (double)image*path.L;
          gauss_sum_old += path.FastExp(-0.5*dist_old*dist_old/sigma2);
          gauss_sum_new += path.FastExp(-0.5*dist_new*dist_new/sigma2);
        }
        gauss_prod_old *= gauss_sum_old;
        gauss_prod_new *= gauss_sum_new;
      }
      old_log_sample_prob += prefactor_sample_prob + log(gauss_prod_old);
      new_log_sample_prob += prefactor_sample_prob + log(gauss_prod_new);

      beadA = beadC;
    }

    // Calculate action change
    double old_action = 0.;
    double new_action = 0.;
    for (auto& action: action_list) {
      // Old action
      path.SetMode(0);
      old_action += action->GetAction(bead0, bead1, particles, level_i);

      // New action
      path.SetMode(1);
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
