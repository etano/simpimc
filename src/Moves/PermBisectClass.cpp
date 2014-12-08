#include "PermBisectClass.h"

void PermBisect::Init(Input &in)
{
  nLevel = in.getAttribute<int>("nLevel");
  int maxPossibleLevel = floor(log2(path.nBead));
  if (nLevel > maxPossibleLevel)
    cout << "Warning: nLevel > maxPossibleLevel!" << endl;
  if (path.PBC)
    nImages = in.getAttribute<int>("nImages");
  else
    nImages = 0;
  species = in.getAttribute<string>("species");
  nPermPart = in.getAttribute<int>("nPermPart");
  epsilon = in.getAttribute<double>("epsilon",1.e-100);
  logEpsilon = log(epsilon);

  // Set species things
  path.GetSpeciesInfo(species,iSpecies);
  nPart = path.speciesList[iSpecies]->nPart;

  // Initialize constant cofactor
  nBisectBeads = 1<<nLevel; // Number of beads in bisection
  lambda = path.speciesList[iSpecies]->lambda;
  i4LambdaTauNBisectBeads = 1./(4.*lambda*path.tau*nBisectBeads);

  // Initiate permutation table
  t.zeros(nPart,nPart);

  // Generate all cycles
  BuildCycles();

  // Initiate acceptance ratio counters
  permAttempt.zeros(nPermType);
  permAccept.zeros(nPermType);
}

// Construct all cycles of particle exchanges
void PermBisect::BuildCycles()
{
  // Generate possible cycles of nPermPart particles
  vector<int> tmpPossCycle;
  for (int i=0; i<nPermPart; ++i)
    tmpPossCycle.push_back(i);
  vector< vector<int> > tmpPossCycles;
  genPerm(tmpPossCycles, tmpPossCycle);
  nPermType = tmpPossCycles.size();
  field<Cycle> possible_cycles;
  possible_cycles.set_size(nPermType);
  for (int i=0; i<nPermType; ++i) {
    mat<int> iP(nPermPart,nPermPart);
    iP.zeros();
    possible_cycles(i).perm.set_size(nPermPart);
    for (int j=0; j<nPermPart; ++j) {
      possible_cycles(i).perm(j) = tmpPossCycles[i][j];
      for (int k=0; k<nPermPart; ++k)
        if (tmpPossCycles[i][k] == j)
          iP(j,k) = 1;
    }
    vec<int> ic(nPermPart);
    for (int j=0; j<nPermPart; ++j)
      ic(j) = tmpPossCycle[j];
    ic = iP * ic;
    possible_cycles(i).iPerm = ic;
  }

  // Run through permatation types
  vector<int> tmpCycle;
  for (int i=0; i<nPart; ++i)
    tmpCycle.push_back(i);
  vector< vector<int> > tmpCycles;
  genCombPermK(tmpCycles, tmpCycle, nPermPart, false, false);
  int nCycle = tmpCycles.size();
  all_cycles.set_size(nPermType * nCycle);
  int permIndex = 0;
  for (unsigned int iPermType=0; iPermType<nPermType; iPermType++) {
    for (int iCycle=0; iCycle<nCycle; ++iCycle) {
      Cycle& c = all_cycles(permIndex);
      c.type = iPermType;
      c.index = permIndex;
      c.part.set_size(tmpCycles.size());
      for (int i=0; i<tmpCycles.size(); ++i)
        c.part(i) = tmpCycles[iCycle][i];
      c.perm = possible_cycles(iPermType).perm;
      c.iPerm = possible_cycles(iPermType).iPerm;
      permIndex += 1;
    }
  }

}

// Accept current move
void PermBisect::Accept()
{
  nAttempt++;
  nAccept++;

  // Accept move, so store things
  for (unsigned int iP=0; iP<nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
    path(iSpecies,iP,bead1) -> storePrev();
    path(iSpecies,iP,bead1-1) -> storeNext();
  }
  assignParticleLabels();
  path.storeR(affBeads);
  path.storeRhoKP(affBeads);
  for (int iB=bead0; iB<=bead1; ++iB)
    path.storeRhoK(iB,iSpecies);

  // Increment permutation counter
  permAttempt(permType) += 1;
  permAccept(permType) += 1;

  // Call accept for each action
  for (auto& action: actionList)
    action->Accept();
}

// Reject current move
void PermBisect::Reject()
{
  nAttempt++;

  // Restore things
  for (unsigned int iP=0; iP<nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
    path(iSpecies,iP,bead1) -> restorePrev();
    path(iSpecies,iP,bead1-1) -> restoreNext();
  }
  path.restoreR(affBeads);
  path.restoreRhoKP(affBeads);
  for (int iB=bead0; iB<=bead1; ++iB)
    path.restoreRhoK(iB,iSpecies);

  // Increment permutation counter
  permAttempt(permType) += 1;

  // Call reject for each action
  for (auto& action: actionList)
    action->Reject();
}

// Perform the permuting bisection
int PermBisect::Attempt()
{
  bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle

  // Set up permutation
  cycles.clear();
  double permTot0 = constructPermTable(); // Permutation weight table
  int cycleIndex = selectCycle(permTot0);
  Cycle* c = cycles[cycleIndex];

  // Set up pointers
  vector< pair<int,int> > particles;
  int nPartPerm = c->part.size();
  field<Bead*> beadI(nPartPerm), beadFm1(nPartPerm), beadF(nPartPerm);
  for (unsigned int i=0; i<nPartPerm; i++) {
    beadI(i) = path(iSpecies,c->part(i),bead0);
    beadFm1(i) = beadI(i)->nextB(nBisectBeads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(std::make_pair(iSpecies,c->part(i)));
  }

  // Permute particles
  permuteBeads(beadFm1, beadF, c);

  // Note affected beads
  field<Bead*> beadA(nPartPerm);
  affBeads.clear();
  for (unsigned int i=0; i<nPartPerm; i++) {
    for(beadA(i) = beadI(i); beadA(i) != beadF(i); beadA(i) = beadA(i) -> next)
      affBeads.push_back(beadA(i));
  }

  // Perform the bisection (move exactly through kinetic action)
  field<Bead*> beadB(nPartPerm), beadC(nPartPerm);
  double oldCycleWeight = -log(c->weight);
  double prevActionChange = oldCycleWeight;
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
    for (unsigned int i=0; i<nPartPerm; i++) {
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

        beadA(i) = beadC(i);
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

  // Have to check this if neighborhood has changed
  if (nPermPart < nPart) {
    // Construct Permutation Table
    path.SetMode(1);
    double permTot1 = constructPermTable();
  
    // Decide whether or not to accept whole bisection
    if ((oldCycleWeight*permTot0/permTot1) < rng.unifRand())
      return 0;
  }

  return 1;
}

// Construct permutation table and probabilities
double PermBisect::constructPermTable()
{
  // Update t
  updatePermTable();

  // Run through permatation types
  double totalWeight = 0.;
  for (unsigned int permIndex=0; permIndex<all_cycles.size(); permIndex++) {
    Cycle& c = all_cycles(permIndex);
    c.weight = 1.;
    for (unsigned int iP=0; iP<c.part.size(); iP++)
      c.weight *= t(c.part(iP),c.part(c.perm(iP)));
    if (c.weight > epsilon) {
      totalWeight += c.weight;
      c.contribution = totalWeight;
      cycles.push_back(&c);
    }
  }
  return totalWeight;
}

void PermBisect::updatePermTable()
{
  // Set initial and final beads
  field<Bead*> b0(nPart), b1(nPart);
  for (unsigned int iP=0; iP<nPart; iP++) {
    b0(iP) = path(iSpecies,iP,bead0);
    b1(iP) = path.GetNextBead(b0(iP),nBisectBeads);
  }

  // Construct t table
  vec<double> dr_ij(path.nD), dr_ii(path.nD);
  double exponent;
  for (unsigned int i=0; i<nPart; i++) {
    path.Dr(b0(i), b1(i), dr_ii);
    for (unsigned int j=0; j<nPart; j++) {
      path.Dr(b0(i), b1(j), dr_ij);
      exponent = (-dot(dr_ij, dr_ij) + dot(dr_ii, dr_ii))*i4LambdaTauNBisectBeads;
      if (exponent > logEpsilon)
        t(i,j) = path.fexp(exponent);
      else
        t(i,j) = 0.;
    }
  }

}

int PermBisect::selectCycle(double permTot)
{
  double x = rng.unifRand(0.,permTot);
  int hi = cycles.size();
  int lo = 0;
  if (x < cycles[0]->contribution)
    return 0;
  while (hi - lo > 1) {
    int mid = (hi+lo)>>1;
    if (x < cycles[mid]->contribution)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

// Permute paths between b0 and b1 given cycle
void PermBisect::permuteBeads(field<Bead*>& b0, field<Bead*>& b1, Cycle* c)
{
  // Set permutation type
  permType = c->type;

  // Execute the permutation
  int nPerm = c->part.size();
  for (unsigned int i=0; i<nPerm; i++)
    b0(i)->next = b1(c->perm(i));
  for (unsigned int i=0; i<nPerm; i++)
    b1(i)->prev = b0(c->iPerm(i));
  for (unsigned int i=0; i<nPerm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisect::assignParticleLabels()
{
  Bead *b;
//  for (unsigned int iP=0; iP<nPart; iP++) {
//    b = path(iSpecies,iP,0);
//    for (unsigned int iB=0; iB<path.nBead; iB++) {
//      path(iSpecies,iP,iB) = b;
//      path(iSpecies,iP,iB)->p = iP;
//      b = b->next;
//    }
//  }
  for (unsigned int iP=0; iP<nPart; iP++) {
    b = path(iSpecies,iP,bead1-1);
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

void PermBisect::Write()
{
  // Write
  if (firstTime) {
    out.Write("/Moves/"+name+"/nPermType", nPermType);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "permAttempt", permAttempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "permAccept", permAccept);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "permAttempt", permAttempt);
    out.AppendDataSet("/Moves/"+name+"/", "permAccept", permAccept);
  }

  // Reset
  permAttempt.zeros();
  permAccept.zeros();

  Move::Write();

}
