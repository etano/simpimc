#include "perm_bisect_class.h"

void PermBisect::Init(Input &in)
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
  n_perm_part = in.GetAttribute<uint>("n_perm_part");
  epsilon = in.GetAttribute<double>("epsilon",1.e-100);
  log_epsilon = log(epsilon);

  // Set species things
  path.GetSpeciesInfo(species,species_i);
  n_part = path.species_list[species_i]->n_part;

  // Generate action list
  std::vector<std::string> species_list;
  species_list.push_back(species);
  GenerateActionList(species_list);

  // Initialize constant cofactor
  n_bisect_beads = 1<<n_level; // Number of beads in bisection
  lambda = path.species_list[species_i]->lambda;
  i4_lambda_tau_n_bisect_beads = 1./(4.*lambda*path.tau*n_bisect_beads);

  // Initiate permutation table
  t.zeros(n_part,n_part);

  // Generate all cycles
  BuildCycles();

  // Initiate acceptance ratio counters
  perm_attempt.zeros(n_perm_type);
  perm_accept.zeros(n_perm_type);
}

// Construct all cycles of particle exchanges
void PermBisect::BuildCycles()
{
  // Generate possible cycles of n_perm_part particles
  std::vector<uint> tmp_poss_cycle;
  for (uint i=0; i<n_perm_part; ++i)
    tmp_poss_cycle.push_back(i);
  std::vector< std::vector<uint>> tmp_poss_cycles;
  GenPerm(tmp_poss_cycles, tmp_poss_cycle);
  n_perm_type = tmp_poss_cycles.size();
  field<Cycle> possible_cycles;
  possible_cycles.set_size(n_perm_type);
  for (uint i=0; i<n_perm_type; ++i) {
    mat<uint> p_i(n_perm_part,n_perm_part);
    p_i.zeros();
    possible_cycles(i).perm.set_size(n_perm_part);
    for (uint j=0; j<n_perm_part; ++j) {
      possible_cycles(i).perm(j) = tmp_poss_cycles[i][j];
      for (uint k=0; k<n_perm_part; ++k)
        if (tmp_poss_cycles[i][k] == j)
          p_i(j,k) = 1;
    }
    vec<uint> ic(n_perm_part);
    for (uint j=0; j<n_perm_part; ++j)
      ic(j) = tmp_poss_cycle[j];
    ic = p_i * ic;
    possible_cycles(i).i_perm = ic;
  }

  // Run through permatation types
  std::vector<uint> tmp_cycle;
  for (uint i=0; i<n_part; ++i)
    tmp_cycle.push_back(i);
  std::vector<std::vector<uint>> tmp_cycles;
  GenCombPermK(tmp_cycles, tmp_cycle, n_perm_part, false, false);
  uint n_cycle = tmp_cycles.size();
  all_cycles.set_size(n_perm_type * n_cycle);
  uint perm_index = 0;
  for (uint perm_type_i=0; perm_type_i<n_perm_type; perm_type_i++) {
    for (uint cycle_i=0; cycle_i<n_cycle; ++cycle_i) {
      Cycle &c = all_cycles(perm_index);
      c.type = perm_type_i;
      c.index = perm_index;
      c.part.set_size(tmp_cycles[cycle_i].size());
      for (uint i=0; i<tmp_cycles[cycle_i].size(); ++i)
        c.part(i) = tmp_cycles[cycle_i][i];
      c.perm = possible_cycles(perm_type_i).perm;
      c.i_perm = possible_cycles(perm_type_i).i_perm;
      perm_index += 1;
    }
  }

}

// Accept current move
void PermBisect::Accept()
{
  n_attempt++;
  n_accept++;

  // Accept move, so store things
  for (uint p_i=0; p_i<n_part; p_i++) { // todo: can make this more efficient by only restoring touched particles
    path(species_i,p_i,bead1)->StorePrev();
    path(species_i,p_i,bead1-1)->StoreNext();
  }
  AssignParticleLabels();
  path.StoreR(affected_beads);
  path.StoreRhoKP(affected_beads);
  for (uint b_i=bead0; b_i<=bead1; ++b_i)
    path.StoreRhoK(b_i,species_i);

  // Increment permutation counter
  perm_attempt(perm_type) += 1;
  perm_accept(perm_type) += 1;

  // Call accept for each action
  for (auto& action: action_list)
    action->Accept();
}

// Reject current move
void PermBisect::Reject()
{
  n_attempt++;

  // Restore things
  for (uint p_i=0; p_i<n_part; p_i++) { // todo: can make this more efficient by only restoring touched particles
    path(species_i,p_i,bead1)->RestorePrev();
    path(species_i,p_i,bead1-1)->RestoreNext();
  }
  path.RestoreR(affected_beads);
  path.RestoreRhoKP(affected_beads);
  for (uint b_i=bead0; b_i<=bead1; ++b_i)
    path.RestoreRhoK(b_i,species_i);

  // Increment permutation counter
  perm_attempt(perm_type) += 1;

  // Call reject for each action
  for (auto& action: action_list)
    action->Reject();
}

// Perform the permuting bisection
bool PermBisect::Attempt()
{
  bead0 = rng.UnifRand(path.n_bead) - 1;  // Pick first bead at random
  bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
  bool roll_over = bead1 > (path.n_bead-1);  // See if bisection overflows to next particle
  // Set up permutation
  cycles.clear();
  double perm_tot_0 = ConstructPermTable(); // Permutation weight table
  uint cycle_index = SelectCycle(perm_tot_0);
  Cycle* c = cycles[cycle_index]; // TODO: Make a smart pointer

  // Set up pointers
  std::vector<std::pair<uint,uint>> particles;
  uint n_perm_part = c->part.size();
  field< std::shared_ptr<Bead>> beadI(n_perm_part), beadFm1(n_perm_part), beadF(n_perm_part);
  for (uint i=0; i<n_perm_part; i++) {
    beadI(i) = path(species_i,c->part(i),bead0);
    beadFm1(i) = beadI(i)->NextB(n_bisect_beads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(std::make_pair(species_i,c->part(i)));
  }

  // Permute particles
  PermuteBeads(beadFm1, beadF, c);

  // Note affected beads
  field< std::shared_ptr<Bead>> beadA(n_perm_part);
  affected_beads.clear();
  for (uint i=0; i<n_perm_part; i++) {
    for(beadA(i) = beadI(i); beadA(i) != beadF(i); beadA(i) = beadA(i) -> next)
      affected_beads.push_back(beadA(i));
  }

  // Perform the bisection (move exactly through kinetic action)
  field< std::shared_ptr<Bead>> beadB(n_perm_part), beadC(n_perm_part);
  double oldCycleWeight = -log(c->weight);
  double prev_action_change = oldCycleWeight;
  double prefactor_sample_prob = 0.;
  vec<double> r_bar_old(path.n_d), delta_old(path.n_d), r_bar_new(path.n_d), delta_new(path.n_d);
  double gauss_prod_old, gauss_sum_old, dist_old, gauss_prod_new, gauss_sum_new, dist_new;
  for (int level_i = n_level-1; level_i >= 0; level_i -= 1) {

    // Level specific quantities
    uint skip = 1<<level_i;
    double levelTau = path.tau*skip;
    double sigma2 = lambda*levelTau;
    double sigma = sqrt(sigma2);

    // Calculate sampling probability
    double old_log_sample_prob = 0.;
    double new_log_sample_prob = 0.;
    for (uint i=0; i<n_perm_part; i++) {
      beadA(i) = beadI(i);
      while(beadA(i) != beadF(i)) {
        // Old sampling
        path.SetMode(0);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        vec<double> r_bar_old(path.RBar(beadC(i), beadA(i)));
        vec<double> delta_old(path.Dr(beadB(i), r_bar_old));

        // New sampling
        path.SetMode(1);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        vec<double> r_bar_new(path.RBar(beadC(i), beadA(i)));
        vec<double> delta_new(path.n_d);
        rng.NormRand(delta_new, 0., sigma);
        path.PutInBox(delta_new);
        beadB(i)->r = r_bar_new + delta_new;

        // Get sampling probs
        gauss_prod_old = 1.;
        gauss_prod_new = 1.;
        for (uint iD=0; iD<path.n_d; iD++) {
          gauss_sum_old = 0.;
          gauss_sum_new = 0.;
          for (int image=-n_images; image<=n_images; image++) {
            dist_old = delta_old(iD) + (double)image*path.L;
            dist_new = delta_new(iD) + (double)image*path.L;
            gauss_sum_old += path.FastExp(-0.5*dist_old*dist_old/sigma2);
            gauss_sum_new += path.FastExp(-0.5*dist_new*dist_new/sigma2);
          }
          gauss_prod_old *= gauss_sum_old;
          gauss_prod_new *= gauss_sum_new;
        }
        old_log_sample_prob += prefactor_sample_prob + log(gauss_prod_old);
        new_log_sample_prob += prefactor_sample_prob + log(gauss_prod_new);

        beadA(i) = beadC(i);
      }
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

    // Calculate acceptance ratio
    double log_sample_ratio = -new_log_sample_prob + old_log_sample_prob;
    double current_action_change = new_action - old_action;
    double log_accept_probablity = log_sample_ratio - current_action_change + prev_action_change;

    // Metropolis step
    if (log_accept_probablity < log(rng.UnifRand()))
      return 0;

    prev_action_change = current_action_change;
  }

  // Have to check this if neighborhood has changed
  if (n_perm_part < n_part) {
    // Construct Permutation Table
    path.SetMode(1);
    double perm_tot_1 = ConstructPermTable();
  
    // Decide whether or not to accept whole bisection
    if ((perm_tot_0/perm_tot_1) < rng.UnifRand())
      return 0;
  }

  return 1;
}

// Construct permutation table and probabilities
double PermBisect::ConstructPermTable()
{
  // Update t
  UpdatePermTable();

  // Run through permatation types
  double total_weight = 0.;
  for (uint perm_index=0; perm_index<all_cycles.size(); perm_index++) {
    Cycle& c = all_cycles(perm_index);
    c.weight = 1.;
    for (uint p_i=0; p_i<c.part.size(); p_i++)
      c.weight *= t(c.part(p_i),c.part(c.perm(p_i)));
    if (c.weight > epsilon) {
      total_weight += c.weight;
      c.contribution = total_weight;
      cycles.push_back(&c);
    }
  }
  return total_weight;
}

void PermBisect::UpdatePermTable()
{
  // Set initial and final beads
  field< std::shared_ptr<Bead>> b0(n_part), b1(n_part);
  for (uint p_i=0; p_i<n_part; p_i++) {
    b0(p_i) = path(species_i,p_i,bead0);
    b1(p_i) = path.GetNextBead(b0(p_i),n_bisect_beads);
  }

  // Construct t table
  for (uint i=0; i<n_part; i++) {
    vec<double> dr_ii(path.Dr(b0(i), b1(i)));
    for (uint j=0; j<n_part; j++) {
      vec<double> dr_ij(path.Dr(b0(i), b1(j)));
      double exponent = (-dot(dr_ij, dr_ij) + dot(dr_ii, dr_ii))*i4_lambda_tau_n_bisect_beads;
      if (exponent > log_epsilon)
        t(i,j) = path.FastExp(exponent);
      else
        t(i,j) = 0.;
    }
  }

}

uint PermBisect::SelectCycle(const double permTot)
{
  double x = rng.UnifRand(0.,permTot);
  uint hi = cycles.size();
  uint lo = 0;
  if (x < cycles[0]->contribution)
    return 0;
  while (hi - lo > 1) {
    uint mid = (hi+lo)>>1;
    if (x < cycles[mid]->contribution)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

// Permute paths between b0 and b1 given cycle
void PermBisect::PermuteBeads(field<std::shared_ptr<Bead>>& b0, field<std::shared_ptr<Bead>>& b1, const Cycle* const c)
{
  // Set permutation type
  perm_type = c->type;

  // Execute the permutation
  uint n_perm = c->part.size();
  for (uint i=0; i<n_perm; i++)
    b0(i)->next = b1(c->perm(i));
  for (uint i=0; i<n_perm; i++)
    b1(i)->prev = b0(c->i_perm(i));
  for (uint i=0; i<n_perm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisect::AssignParticleLabels()
{
//  for (uint p_i=0; p_i<n_part; p_i++) {
//    b = path(species_i,p_i,0);
//    for (uint b_i=0; b_i<path.n_bead; b_i++) {
//      path(species_i,p_i,b_i) = b;
//      path(species_i,p_i,b_i)->p = p_i;
//      b = b->next;
//    }
//  }
  for (uint p_i=0; p_i<n_part; p_i++) {
    std::shared_ptr<Bead> b(path(species_i,p_i,bead1-1));
    for (uint b_i=path.bead_loop(bead1-1); b_i<path.n_bead; b_i++) {
      path.species_list[species_i]->bead(p_i,b_i) = b;
      path(species_i,p_i,b_i)->p = p_i;
      b = b->next;
    }
  }

  //for (uint p_i=0; p_i<n_part; p_i++) {
  //  for (uint b_i=0; b_i<path.n_bead; b_i++) {
  //    cout << p_i << " " << b_i << "   " << path(species_i,p_i,b_i)->prev->p << " " << path(species_i,p_i,b_i)->p << " " << path(species_i,p_i,b_i)->next->p << "   " << path(species_i,p_i,b_i)->prev->b << " " << path(species_i,p_i,b_i)->b << " " << path(species_i,p_i,b_i)->next->b << endl;
  //  }
  //}
}

void PermBisect::Write()
{
  // Write
  if (first_time) {
    out.Write("/Moves/"+name+"/n_perm_type", n_perm_type);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
    out.AppendDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
  }

  // Reset
  perm_attempt.zeros();
  perm_accept.zeros();

  Move::Write();

}
