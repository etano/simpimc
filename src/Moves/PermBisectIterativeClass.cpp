#include "PermBisectIterativeClass.h"

void PermBisectIterative::Init(Input &in)
{
  // Read in things
  nLevel = in.getAttribute<int>("nLevel");
  int maxPossibleLevel = floor(log2(path.nBead));
  if (nLevel > maxPossibleLevel)
    cout << "Warning: nLevel > maxPossibleLevel!" << endl;
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;
  species = in.getAttribute<string>("species");
  epsilon = in.getAttribute<double>("epsilon",1.e-100);
  logEpsilon = log(epsilon);

  // Adaptive bisection level
  adaptive = in.getAttribute<int>("adaptive",0);
  if (adaptive)
    targetRatio = in.getAttribute<double>("targetRatio");

  // Set species things
  path.GetSpeciesInfo(species,iSpecies);
  nPart = path.speciesList[iSpecies]->nPart;

  // Generate action list
  std::vector<std::string> speciesList;
  speciesList.push_back(species);
  GenerateActionList(speciesList);

  // Initialize constant cofactor
  nBisectBeads = 1<<nLevel; // Number of beads in bisection
  lambda = path.speciesList[iSpecies]->lambda;
  i4LambdaTauNBisectBeads = 1./(4.*lambda*path.tau*nBisectBeads);

  // Initiate permutation table
  t.zeros(nPart,nPart);

  // Initiate acceptance ratio counters
  permAttempt.zeros(nPart);
  permAccept.zeros(nPart);
  refAccept = 0;
  refAttempt = 0;
}

void PermBisectIterative::Accept()
{
  permAttempt(nPermPart-1) += 1;
  permAccept(nPermPart-1) += 1;

  // Change sign weight for fermions
  if (!(nPermPart%2) && path.speciesList[iSpecies]->fermi)
    path.sign *= -1;

  // Accept move, so store things
  for (unsigned int iP=0; iP<nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
    path(iSpecies,iP,bead1) -> storePrev();
    path(iSpecies,iP,bead1-1) -> storeNext();
  }
  if (nPermPart > 1) // only need to reassign particle labels if actual permutation
    assignParticleLabels();
  path.storeR(affBeads);
  path.storeRhoKP(affBeads);
  for (int iB=bead0; iB<bead1; ++iB)
    path.storeRhoK(iB,iSpecies);

  // Call reject for each action
  for (auto& action: actionList)
    action->Accept();
}

void PermBisectIterative::Reject()
{
  // No need to do some things if bisection isn't attempted
  if (nPermPart > 0) {
    permAttempt(nPermPart-1) += 1;

    // Restore things
    for (unsigned int iP=0; iP<nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
      path(iSpecies,iP,bead1) -> restorePrev();
      path(iSpecies,iP,bead1-1) -> restoreNext();
    }
    path.restoreR(affBeads);
    path.restoreRhoKP(affBeads);
    for (int iB=bead0; iB<bead1; ++iB)
      path.restoreRhoK(iB,iSpecies);
  }

  // Call reject for each action
  for (auto& action: actionList)
    action->Reject();
}

void PermBisectIterative::Reset()
{
  if (adaptive) {
    double acceptRatio = (double) nAccept / (double) nAttempt;
    if (acceptRatio < targetRatio && nLevel > 1)
      nLevel--;
    else if (1<<nLevel <= path.nBead/2)
      nLevel++;
    nBisectBeads = 1<<nLevel; // Number of beads in bisection
    lambda = path.speciesList[iSpecies]->lambda;
    i4LambdaTauNBisectBeads = 1./(4.*lambda*path.tau*nBisectBeads);
  }

  refAccept = 0;
  refAttempt = 0;

  Move::Reset();
}

// Perform the permuting bisection
int PermBisectIterative::Attempt()
{
  bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle
  bool includeRef = path.speciesList[iSpecies]->fixedNode &&
                    ((bead0<=path.refBead && bead1>=path.refBead) ||
                    (rollOver && path.beadLoop[bead1]>=path.refBead));
  if (includeRef)
    refAttempt++;

  // Set up permutation
  Cycle c;
  nPermPart = 0; // reset to indicate if bisection is atempted or not
  if (!selectCycleIterative(c))
    return 0; // do not attempt bisection since permutation not accepted
  nPermPart = c.part.size();

  // Set up pointers
  vector< pair<int,int> > particles;
  field< std::shared_ptr<Bead> > beadI(nPermPart), beadFm1(nPermPart), beadF(nPermPart);
  for (unsigned int i=0; i<nPermPart; i++) {
    beadI(i) = path(iSpecies,c.part(i),bead0);
    beadFm1(i) = beadI(i)->nextB(nBisectBeads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(std::make_pair(iSpecies,c.part(i)));
  }
  if (rollOver) {
    for (int i=0; i<nPermPart; ++i)
      particles.push_back(std::make_pair(iSpecies,beadF(i)->p));
    sort(particles.begin(), particles.end());
    particles.erase(unique(particles.begin(), particles.end()), particles.end());
  }

  // Permute particles
  permuteBeads(beadFm1, beadF, c);

  // Note affected beads
  field< std::shared_ptr<Bead> > beadA(nPermPart);
  affBeads.clear();
  for (unsigned int i=0; i<nPermPart; i++) {
    for(beadA(i) = beadI(i)->next; beadA(i) != beadF(i); beadA(i) = beadA(i)->next)
      affBeads.push_back(beadA(i));
  }

  // Perform the bisection (move exactly through kinetic action)
  field< std::shared_ptr<Bead> > beadB(nPermPart), beadC(nPermPart);
  double prevActionChange = -log(c.weight);
  double prefactorOfSampleProb = 0.;
  vec<double> rBarOld(path.nD), deltaOld(path.nD), rBarNew(path.nD), deltaNew(path.nD);
  double gaussProdOld, gaussSumOld, distOld, gaussProdNew, gaussSumNew, distNew;
  for (int iLevel = nLevel-1; iLevel >= 0; iLevel -= 1) {

    // Level specific quantities
    int skip = 1<<iLevel;
    double levelTau = path.tau*skip;
    double sigma2 = lambda*levelTau;
    double sigma = sqrt(sigma2);

    // Calculate sampling probability
    double oldLogSampleProb = 0.;
    double newLogSampleProb = 0.;
    for (unsigned int i=0; i<nPermPart; i++) {
      beadA(i) = beadI(i);
      while(beadA(i) != beadF(i)) {
        // Old sampling
        path.SetMode(0);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        path.RBar(beadC(i), beadA(i), rBarOld);
        path.Dr(beadB(i), rBarOld, deltaOld);

        // New sampling
        path.SetMode(1);
        beadB(i) = path.GetNextBead(beadA(i),skip);
        beadC(i) = path.GetNextBead(beadB(i),skip);
        path.RBar(beadC(i), beadA(i), rBarNew);
        rng.normRand(deltaNew, 0., sigma);
        path.PutInBox(deltaNew);
        beadB(i)->r = rBarNew + deltaNew;

        // Get sampling probs
        gaussProdOld = 1.;
        gaussProdNew = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          gaussSumOld = 0.;
          gaussSumNew = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            distOld = deltaOld(iD) + (double)image*path.L;
            distNew = deltaNew(iD) + (double)image*path.L;
            gaussSumOld += path.fexp(-0.5*distOld*distOld/sigma2);
            gaussSumNew += path.fexp(-0.5*distNew*distNew/sigma2);
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
    double oldAction = 0.;
    double newAction = 0.;
    for (auto& action: actionList) {
      // Old action
      path.SetMode(0);
      oldAction += action->GetAction(bead0, bead1, particles, iLevel);

      // New action
      path.SetMode(1);
      newAction += action->GetAction(bead0, bead1, particles, iLevel);
    }

    // Calculate acceptance ratio
    double logSampleRatio = -newLogSampleProb + oldLogSampleProb;
    double currActionChange = newAction - oldAction;
    double logAcceptProb = logSampleRatio - currActionChange + prevActionChange;

    // Metropolis step
    if (logAcceptProb < log(rng.unifRand()))
      return 0;

    prevActionChange = currActionChange;
  }

  if (includeRef)
    refAccept++;

  return 1;
}

void PermBisectIterative::updatePermTable()
{
  // Set initial and final beads
  field< std::shared_ptr<Bead> > b0(nPart), b1(nPart);
  for (unsigned int iP=0; iP<nPart; iP++) {
    b0(iP) = path(iSpecies,iP,bead0);
    b1(iP) = b0(iP) -> nextB(nBisectBeads);
  }

  // Construct t table
  vec<double> dr_ij(path.nD), dr_ii(path.nD);
  double exponent;
  for (unsigned int i=0; i<nPart; i++) {
    for (unsigned int j=0; j<nPart; j++) {
      path.Dr(b0(i), b1(j), dr_ij);
      exponent = (-dot(dr_ij,dr_ij))*i4LambdaTauNBisectBeads;
      if (exponent > logEpsilon)
        t(i,j) = path.fexp(exponent);
      else
        t(i,j) = 0.;
    }
  }

}

int PermBisectIterative::selectCycleIterative(Cycle& c)
{
  // Update t
  updatePermTable();
  mat<double> t_c = t;

  // Choose particles
  int p0 = rng.unifRand(nPart) - 1;  // Pick first particle at random
  int p = p0;
  vector<int> ps;
  do {
    // Add particle to ps
    ps.push_back(p);

    // Make sure returning to previous particles is not an option
    for (int i=0; i<ps.size(); ++i)
      t_c(p,ps[i]) = 0.;

    // Allow cycle to close
    if (ps.size() > 0)
      t_c(p,p0) = t(p,p0);

    // Disallow even permutations for fixed-node fermions
    if (path.speciesList[iSpecies]->fermi && path.speciesList[iSpecies]->fixedNode && !(ps.size()%2))
      t_c(p,p0) = 0.;

    // Calculate row total
    double Q_p = 0.;
    double Q_p_c = 0.;
    for (int i=0; i<nPart; ++i) {
      Q_p += t(p,i);
      Q_p_c += t_c(p,i);
    }

    // Decide whether or not to continue
    if ((Q_p_c/Q_p) < rng.unifRand())
      return 0;

    // Select next particle with bisective search
    double x = rng.unifRand();
    double t_Q = 0.;
    for (int i=0; i<nPart; ++i) { // fixme: not doing bisection
      t_Q += t_c(p,i)/Q_p_c;
      if (t_Q > x) {
        p = i;
        break;
      }
    }

  } while (p != p0);

  // Set weight
  c.weight = 1.;
  int nPerm = ps.size();
  for (int i=0; i<nPerm-1; ++i)
    c.weight *= t(ps[i],ps[i+1])/t(ps[i],ps[i]);
  c.weight *= t(ps[nPerm-1],ps[0])/t(ps[nPerm-1],ps[nPerm-1]);

  // Set particles
  c.part.set_size(nPerm);
  for (int i=0; i<nPerm; ++i)
    c.part(i) = ps[i];

  // Set perms
  c.perm.set_size(nPerm);
  for (int i=0; i<nPerm-1; ++i)
    c.perm(i) = i+1;
  c.perm(nPerm-1) = 0;
  c.iPerm.set_size(nPerm);
  c.iPerm(0) = nPerm-1;
  for (int i=1; i<nPerm; ++i)
    c.iPerm(i) = i-1;


  return 1;
}

// Permute paths between b0 and b1 given cycle
void PermBisectIterative::permuteBeads(field< std::shared_ptr<Bead> > &b0, field< std::shared_ptr<Bead> > &b1, Cycle& c)
{
  // Execute the permutation
  int nPerm = c.part.size();
  for (unsigned int i=0; i<nPerm; i++)
    b0(i)->next = b1(c.perm(i));
  for (unsigned int i=0; i<nPerm; i++)
    b1(i)->prev = b0(c.iPerm(i));
  for (unsigned int i=0; i<nPerm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisectIterative::assignParticleLabels()
{
  for (unsigned int iP=0; iP<nPart; iP++) {
    std::shared_ptr<Bead> b(path(iSpecies,iP,bead1-1));
    for (unsigned int iB=path.beadLoop(bead1-1); iB<path.nBead; iB++) {
      path.speciesList[iSpecies]->bead(iP,iB) = b; // fixme: make cleaner using operator
      path(iSpecies,iP,iB)->p = iP;
      b = b->next;
    }
  }

  //for (unsigned int iP=0; iP<nPart; iP++) {
  //  for (unsigned int iB=0; iB<path.nBead; iB++) {
  //    cout << iP << " " << iB << "   " << path(iSpecies,iP,iB)->prev->p << " " << path(iSpecies,iP,iB)->p << " " << path(iSpecies,iP,iB)->next->p << "   " << path(iSpecies,iP,iB)->prev->b << " " << path(iSpecies,iP,iB)->b << " " << path(iSpecies,iP,iB)->next->b << endl;
  //  }
  //}
}

void PermBisectIterative::Write()
{
  // Write
  if (firstTime) {
    out.CreateExtendableDataSet("/Moves/"+name+"/", "permAttempt", permAttempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "permAccept", permAccept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "refAccept", refAccept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "refAttempt", refAttempt);
    if (adaptive)
      out.CreateExtendableDataSet("/Moves/"+name+"/", "nLevel", nLevel);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "permAttempt", permAttempt);
    out.AppendDataSet("/Moves/"+name+"/", "permAccept", permAccept);
    out.AppendDataSet("/Moves/"+name+"/", "refAttempt", refAttempt);
    out.AppendDataSet("/Moves/"+name+"/", "refAccept", refAccept);
    if (adaptive)
      out.AppendDataSet("/Moves/"+name+"/", "nLevel", nLevel);
  }

  // Reset
  permAttempt.zeros();
  permAccept.zeros();

  Move::Write();

}
