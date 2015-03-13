#include "perm_bisect_iterative_class.h"

void PermBisectIterative::Init(Input &in)
{
  // Read in things
  n_level = in.GetAttribute<uint>("n_level");
  uint max_possible_level = floor(log2(path.n_bead));
  if (n_level > max_possible_level)
    std::cout << "Warning: n_level > max_possible_level!" << std::endl;
  if (path.pbc)
    n_images = in.GetAttribute<int>("n_images");
  else
    n_images = 0;
  species = in.GetAttribute<std::string>("species");
  epsilon = in.GetAttribute<double>("epsilon",1.e-100);
  log_epsilon = log(epsilon);

  // Adaptive bisection level
  adaptive = in.GetAttribute<bool>("adaptive",0);
  if (adaptive)
    target_ratio = in.GetAttribute<double>("target_ratio");

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

  // Initiate acceptance ratio counters
  perm_attempt.zeros(n_part);
  perm_accept.zeros(n_part);
  ref_accept = 0;
  ref_attempt = 0;
}

void PermBisectIterative::Accept()
{
  perm_attempt(n_perm_part-1) += 1;
  perm_accept(n_perm_part-1) += 1;

  // Change sign weight for fermions
  if (!(n_perm_part%2) && path.species_list[species_i]->fermi)
    path.sign *= -1;

  // Accept move, so store things
  for (uint p_i=0; p_i<n_part; p_i++) { // todo: can make this more efficient by only restoring touched particles
    path(species_i,p_i,bead1)->StorePrev();
    path(species_i,p_i,bead1-1)->StoreNext();
  }
  if (n_perm_part > 1) // only need to reassign particle labels if actual permutation
    AssignParticleLabels();
  path.StoreR(affected_beads);
  path.StoreRhoKP(affected_beads);
  for (uint b_i=bead0; b_i<bead1; ++b_i)
    path.StoreRhoK(b_i,species_i);

  // Call reject for each action
  for (auto& action: action_list)
    action->Accept();
}

void PermBisectIterative::Reject()
{
  // No need to do some things if bisection isn't attempted
  if (n_perm_part > 0) {
    perm_attempt(n_perm_part-1) += 1;

    // Restore things
    for (uint p_i=0; p_i<n_part; p_i++) { // todo: can make this more efficient by only restoring touched particles
      path(species_i,p_i,bead1)->RestorePrev();
      path(species_i,p_i,bead1-1)->RestoreNext();
    }
    path.RestoreR(affected_beads);
    path.RestoreRhoKP(affected_beads);
    for (uint b_i=bead0; b_i<bead1; ++b_i)
      path.RestoreRhoK(b_i,species_i);
  }

  // Call reject for each action
  for (auto& action: action_list)
    action->Reject();
}

void PermBisectIterative::Reset()
{
  if (adaptive) {
    double acceptRatio = (double) n_accept / (double) n_attempt;
    if (acceptRatio < target_ratio && n_level > 1)
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

// Perform the permuting bisection
bool PermBisectIterative::Attempt()
{
  bead0 = rng.UnifRand(path.n_bead) - 1;  // Pick first bead at random
  bead1 = bead0 + n_bisect_beads; // Set last bead in bisection
  roll_over = bead1 > (path.n_bead-1);  // See if bisection overflows to next particle
  bool includeRef = path.species_list[species_i]->fixed_node &&
                    ((bead0<=path.ref_bead && bead1>=path.ref_bead) ||
                    (roll_over && path.bead_loop[bead1]>=path.ref_bead));
  if (includeRef)
    ref_attempt++;

  // Set up permutation
  Cycle c;
  n_perm_part = 0; // reset to indicate if bisection is atempted or not
  if (!SelectCycleIterative(c))
    return 0; // do not attempt bisection since permutation not accepted
  n_perm_part = c.part.size();

  // Set up pointers
  std::vector<std::pair<uint,uint>> particles;
  field< std::shared_ptr<Bead>> beadI(n_perm_part), beadFm1(n_perm_part), beadF(n_perm_part);
  for (uint i=0; i<n_perm_part; i++) {
    beadI(i) = path(species_i,c.part(i),bead0);
    beadFm1(i) = beadI(i)->NextB(n_bisect_beads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(std::make_pair(species_i,c.part(i)));
  }
  if (roll_over) {
    for (uint i=0; i<n_perm_part; ++i)
      particles.push_back(std::make_pair(species_i,beadF(i)->p));
    sort(particles.begin(), particles.end());
    particles.erase(unique(particles.begin(), particles.end()), particles.end());
  }

  // Permute particles
  PermuteBeads(beadFm1, beadF, c);

  // Note affected beads
  field< std::shared_ptr<Bead>> beadA(n_perm_part);
  affected_beads.clear();
  for (uint i=0; i<n_perm_part; i++) {
    for(beadA(i) = beadI(i)->next; beadA(i) != beadF(i); beadA(i) = beadA(i)->next)
      affected_beads.push_back(beadA(i));
  }

  // Perform the bisection (move exactly through kinetic action)
  field< std::shared_ptr<Bead>> beadB(n_perm_part), beadC(n_perm_part);
  double prevActionChange = -log(c.weight);
  double prefactorOfSampleProb = 0.;
  vec<double> rBarOld(path.n_d), deltaOld(path.n_d), rBarNew(path.n_d), deltaNew(path.n_d);
  double gaussProdOld, gaussSumOld, distOld, gaussProdNew, gaussSumNew, distNew;
  for (int iLevel = n_level-1; iLevel >= 0; iLevel -= 1) {

    // Level specific quantities
    uint skip = 1<<iLevel;
    double levelTau = path.tau*skip;
    double sigma2 = lambda*levelTau;
    double sigma = sqrt(sigma2);

    // Calculate sampling probability
    double oldLogSampleProb = 0.;
    double newLogSampleProb = 0.;
    for (uint i=0; i<n_perm_part; i++) {
      beadA(i) = beadI(i);
      while(beadA(i) != beadF(i)) {
        // Old sampling
        path.SetMode(0);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        vec<double> rBarOld(path.RBar(beadC(i), beadA(i)));
        vec<double> deltaOld(path.Dr(beadB(i), rBarOld));

        // New sampling
        path.SetMode(1);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        vec<double> rBarNew(path.RBar(beadC(i), beadA(i)));
        vec<double> deltaNew(path.n_d);
        rng.NormRand(deltaNew, 0., sigma);
        path.PutInBox(deltaNew);
        beadB(i)->r = rBarNew + deltaNew;

        // Get sampling probs
        gaussProdOld = 1.;
        gaussProdNew = 1.;
        for (uint d_i=0; d_i<path.n_d; d_i++) {
          gaussSumOld = 0.;
          gaussSumNew = 0.;
          for (int image=-n_images; image<=n_images; image++) {
            distOld = deltaOld(d_i) + (double)image*path.L;
            distNew = deltaNew(d_i) + (double)image*path.L;
            gaussSumOld += path.FastExp(-0.5*distOld*distOld/sigma2);
            gaussSumNew += path.FastExp(-0.5*distNew*distNew/sigma2);
          }
          gaussProdOld *= gaussSumOld;
          gaussProdNew *= gaussSumNew;
        }
        oldLogSampleProb += prefactorOfSampleProb + log(gaussProdOld);
        newLogSampleProb += prefactorOfSampleProb + log(gaussProdNew);

        path.SetMode(1);
        beadA(i) = path.GetNextBead(beadA(i),2*skip);
      }
    }

    // Calculate action change
    double old_action = 0.;
    double new_action = 0.;
    for (auto& action: action_list) {
      // Old action
      path.SetMode(0);
      old_action += action->GetAction(bead0, bead1, particles, iLevel);

      // New action
      path.SetMode(1);
      new_action += action->GetAction(bead0, bead1, particles, iLevel);
    }

    // Calculate acceptance ratio
    double logSampleRatio = -newLogSampleProb + oldLogSampleProb;
    double current_action_change = new_action - old_action;
    double log_accept_probablity = logSampleRatio - current_action_change + prevActionChange;

    // Metropolis step
    if (log_accept_probablity < log(rng.UnifRand()))
      return 0;

    prevActionChange = current_action_change;
  }

  if (includeRef)
    ref_accept++;

  return 1;
}

void PermBisectIterative::UpdatePermTable()
{
  // Set initial and final beads
  field< std::shared_ptr<Bead>> b0(n_part), b1(n_part);
  for (uint p_i=0; p_i<n_part; p_i++) {
    b0(p_i) = path(species_i,p_i,bead0);
    b1(p_i) = b0(p_i)->NextB(n_bisect_beads);
  }

  // Construct t table
  #pragma omp parallel for collapse(2)
  for (uint i=0; i<n_part; i++) {
    for (uint j=i; j<n_part; j++) {
      vec<double> dr_ij(path.Dr(b0(i), b1(j)));
      double exponent = (-dot(dr_ij,dr_ij))*i4_lambda_tau_n_bisect_beads;
      if (exponent > log_epsilon)
        t(i,j) = path.FastExp(exponent);
      else
        t(i,j) = 0.;
      t(j,i) = t(i,j);
    }
  }

}

bool PermBisectIterative::SelectCycleIterative(Cycle& c)
{
  // Update t
  UpdatePermTable();
  mat<double> t_c = t;

  // Choose particles
  uint p0 = rng.UnifRand(n_part) - 1;  // Pick first particle at random
  uint p = p0;
  std::vector<uint> ps;
  do {
    // Add particle to ps
    ps.push_back(p);

    // Make sure returning to previous particles is not an option
    for (uint i=0; i<ps.size(); ++i)
      t_c(p,ps[i]) = 0.;

    // Allow cycle to close
    if (ps.size() > 0)
      t_c(p,p0) = t(p,p0);

    // Disallow even permutations for fixed-node fermions
    if (path.species_list[species_i]->fermi && path.species_list[species_i]->fixed_node && !(ps.size()%2))
      t_c(p,p0) = 0.;

    // Calculate row total
    double Q_p = 0.;
    double Q_p_c = 0.;
    for (uint i=0; i<n_part; ++i) {
      Q_p += t(p,i);
      Q_p_c += t_c(p,i);
    }

    // Decide whether or not to continue
    if ((Q_p_c/Q_p) < rng.UnifRand())
      return 0;

    // Select next particle with bisective search
    double x = rng.UnifRand();
    double t_Q = 0.;
    for (uint i=0; i<n_part; ++i) { // fixme: not doing bisection
      t_Q += t_c(p,i)/Q_p_c;
      if (t_Q > x) {
        p = i;
        break;
      }
    }

  } while (p != p0);

  // Set weight
  c.weight = 1.;
  uint nPerm = ps.size();
  for (uint i=0; i<nPerm-1; ++i)
    c.weight *= t(ps[i],ps[i+1])/t(ps[i],ps[i]);
  c.weight *= t(ps[nPerm-1],ps[0])/t(ps[nPerm-1],ps[nPerm-1]);

  // Set particles
  c.part.set_size(nPerm);
  for (uint i=0; i<nPerm; ++i)
    c.part(i) = ps[i];

  // Set perms
  c.perm.set_size(nPerm);
  for (uint i=0; i<nPerm-1; ++i)
    c.perm(i) = i+1;
  c.perm(nPerm-1) = 0;
  c.i_perm.set_size(nPerm);
  c.i_perm(0) = nPerm-1;
  for (uint i=1; i<nPerm; ++i)
    c.i_perm(i) = i-1;


  return 1;
}

// Permute paths between b0 and b1 given cycle
void PermBisectIterative::PermuteBeads(field<std::shared_ptr<Bead>>& b0, field<std::shared_ptr<Bead>>& b1, const Cycle& c)
{
  // Execute the permutation
  uint nPerm = c.part.size();
  for (uint i=0; i<nPerm; i++)
    b0(i)->next = b1(c.perm(i));
  for (uint i=0; i<nPerm; i++)
    b1(i)->prev = b0(c.i_perm(i));
  for (uint i=0; i<nPerm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisectIterative::AssignParticleLabels()
{
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

void PermBisectIterative::Write()
{
  // Write
  if (first_time) {
    out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
    if (adaptive)
      out.CreateExtendableDataSet("/Moves/"+name+"/", "n_level", n_level);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
    out.AppendDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
    out.AppendDataSet("/Moves/"+name+"/", "ref_attempt", ref_attempt);
    out.AppendDataSet("/Moves/"+name+"/", "ref_accept", ref_accept);
    if (adaptive)
      out.AppendDataSet("/Moves/"+name+"/", "n_level", n_level);
  }

  // Reset
  perm_attempt.zeros();
  perm_accept.zeros();

  Move::Write();

}
