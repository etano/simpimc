#include "perm_bisect_iterative_class.h"

// Perform the permuting bisection
bool PermBisectIterative::Attempt()
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
  n_perm_part = 0; // reset to indicate if bisection is atempted or not
  if (!SelectCycleIterative(c))
    return 0; // do not attempt bisection since permutation not accepted
  n_perm_part = c.part.size();

  // Set up pointers
  std::vector<std::pair<uint32_t,uint32_t>> particles;
  field< std::shared_ptr<Bead>> beadI(n_perm_part), beadFm1(n_perm_part), beadF(n_perm_part);
  for (uint32_t i=0; i<n_perm_part; i++) {
    beadI(i) = path(species_i,c.part(i),bead0);
    beadFm1(i) = beadI(i)->NextB(n_bisect_beads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(std::make_pair(species_i,c.part(i)));
  }
  if (roll_over) {
    for (uint32_t i=0; i<n_perm_part; ++i)
      particles.push_back(std::make_pair(species_i,beadF(i)->p));
    sort(particles.begin(), particles.end());
    particles.erase(unique(particles.begin(), particles.end()), particles.end());
  }

  // Permute particles
  PermuteBeads(beadFm1, beadF, c);

  // Note affected beads
  field< std::shared_ptr<Bead>> beadA(n_perm_part);
  affected_beads.clear();
  for (uint32_t i=0; i<n_perm_part; i++) {
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
    uint32_t skip = 1<<iLevel;
    double levelTau = path.tau*skip;
    double sigma2 = lambda*levelTau;
    double sigma = sqrt(sigma2);

    // Calculate sampling probability
    double oldLogSampleProb = 0.;
    double newLogSampleProb = 0.;
    for (uint32_t i=0; i<n_perm_part; i++) {
      beadA(i) = beadI(i);
      while(beadA(i) != beadF(i)) {
        // Old sampling
        path.SetMode(OLD_MODE);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        vec<double> rBarOld(path.RBar(beadC(i), beadA(i)));
        vec<double> deltaOld(path.Dr(beadB(i), rBarOld));

        // New sampling
        path.SetMode(NEW_MODE);
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
        for (uint32_t d_i=0; d_i<path.n_d; d_i++) {
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

        path.SetMode(NEW_MODE);
        beadA(i) = path.GetNextBead(beadA(i),2*skip);
      }
    }

    // Calculate action change
    double old_action = 0.;
    double new_action = 0.;
    for (auto& action: action_list) {
      // Old action
      path.SetMode(OLD_MODE);
      old_action += action->GetAction(bead0, bead1, particles, iLevel);

      // New action
      path.SetMode(NEW_MODE);
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

  if (include_ref)
    ref_accept++;

  return 1;
}

void PermBisectIterative::UpdatePermTable()
{
  // Set initial and final beads
  field< std::shared_ptr<Bead>> b0(n_part), b1(n_part);
  for (uint32_t p_i=0; p_i<n_part; p_i++) {
    b0(p_i) = path(species_i,p_i,bead0);
    b1(p_i) = b0(p_i)->NextB(n_bisect_beads);
  }

  // Construct t table
  for (uint32_t i=0; i<n_part; i++) {
    for (uint32_t j=i; j<n_part; j++) {
      vec<double> dr_ij(path.Dr(b0(i), b1(j)));
      double exponent = (-dot(dr_ij,dr_ij))*i_4_lambda_tau_n_bisect_beads;
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
  uint32_t p0 = rng.UnifRand(n_part) - 1;  // Pick first particle at random
  uint32_t p = p0;
  std::vector<uint32_t> ps;
  do {
    // Add particle to ps
    ps.push_back(p);

    // Make sure returning to previous particles is not an option
    for (uint32_t i=0; i<ps.size(); ++i)
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
    for (uint32_t i=0; i<n_part; ++i) {
      Q_p += t(p,i);
      Q_p_c += t_c(p,i);
    }

    // Decide whether or not to continue
    if ((Q_p_c/Q_p) < rng.UnifRand())
      return 0;

    // Select next particle with bisective search
    double x = rng.UnifRand();
    double t_Q = 0.;
    for (uint32_t i=0; i<n_part; ++i) { // fixme: not doing bisection
      t_Q += t_c(p,i)/Q_p_c;
      if (t_Q > x) {
        p = i;
        break;
      }
    }

  } while (p != p0);

  // Set weight
  c.weight = 1.;
  uint32_t nPerm = ps.size();
  for (uint32_t i=0; i<nPerm-1; ++i)
    c.weight *= t(ps[i],ps[i+1])/t(ps[i],ps[i]);
  c.weight *= t(ps[nPerm-1],ps[0])/t(ps[nPerm-1],ps[nPerm-1]);

  // Set particles
  c.part.set_size(nPerm);
  for (uint32_t i=0; i<nPerm; ++i)
    c.part(i) = ps[i];
  c.type = nPerm;

  // Set perms
  c.perm.set_size(nPerm);
  for (uint32_t i=0; i<nPerm-1; ++i)
    c.perm(i) = i+1;
  c.perm(nPerm-1) = 0;
  c.i_perm.set_size(nPerm);
  c.i_perm(0) = nPerm-1;
  for (uint32_t i=1; i<nPerm; ++i)
    c.i_perm(i) = i-1;


  return 1;
}
