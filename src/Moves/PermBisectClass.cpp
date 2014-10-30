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
  epsilon = in.getAttribute<RealType>("epsilon",1.e-100);
  logEpsilon = log(epsilon);

  // Set species things
  path.GetSpeciesInfo(species,iSpecies,offset);
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
    Imatrix iP(nPermPart,nPermPart);
    iP.zeros();
    possible_cycles(i).perm = tmpPossCycles[i];
    for (int j=0; j<nPermPart; ++j)
      for (int k=0; k<nPermPart; ++k)
        if (tmpPossCycles[i][k] == j)
          iP(j,k) = 1;
    Ivector ic(tmpPossCycle);
    ic = iP * ic;
    possible_cycles(i).iPerm = ic;
  }

  // Run through permatation types
  vector<int> tmpCycle;
  for (int i=offset; i<offset+nPart; ++i)
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
      c.part = tmpCycles[iCycle];
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
  for (unsigned int iP=offset; iP<offset+nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
    path.bead(iP,path.beadLoop(bead1)) -> storePrev();
    path.bead(iP,path.beadLoop(bead1-1)) -> storeNext();
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
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Accept();
}

// Reject current move
void PermBisect::Reject()
{
  nAttempt++;

  // Restore things
  for (unsigned int iP=offset; iP<offset+nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
    path.bead(iP,path.beadLoop(bead1)) -> restorePrev();
    path.bead(iP,path.beadLoop(bead1-1)) -> restoreNext();
  }
  path.restoreR(affBeads);
  path.restoreRhoKP(affBeads);
  for (int iB=bead0; iB<=bead1; ++iB)
    path.restoreRhoK(iB,iSpecies);

  // Increment permutation counter
  permAttempt(permType) += 1;

  // Call reject for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Reject();
}

// Perform the permuting bisection
int PermBisect::Attempt()
{
  bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle

  // Set up permutation
  cycles.clear();
  RealType permTot0 = constructPermTable(); // Permutation weight table
  int cycleIndex = selectCycle(permTot0);
  Cycle* c = cycles[cycleIndex];

  // Set up pointers
  vector<int> particles;
  int nPartPerm = c->part.size();
  field<Bead*> beadI(nPartPerm), beadFm1(nPartPerm), beadF(nPartPerm);
  for (unsigned int i=0; i<nPartPerm; i++) {
    beadI(i) = path.bead(c->part(i),bead0);
    beadFm1(i) = beadI(i)->nextB(nBisectBeads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(c->part(i));
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
  RealType oldCycleWeight = -log(c->weight);
  RealType prevActionChange = oldCycleWeight;
  RealType prefactorOfSampleProb = 0.;
  Tvector rBarOld(path.nD), deltaOld(path.nD), rBarNew(path.nD), deltaNew(path.nD);
  RealType gaussProdOld, gaussSumOld, distOld, gaussProdNew, gaussSumNew, distNew;
  for (int iLevel = nLevel-1; iLevel >= 0; iLevel -= 1) {

    // Level specific quantities
    int skip = 1<<iLevel;
    RealType levelTau = path.tau*skip;
    RealType sigma2 = lambda*levelTau;
    RealType sigma = sqrt(sigma2);

    // Calculate sampling probability
    RealType oldLogSampleProb = 0.;
    RealType newLogSampleProb = 0.;
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
        rng.normRand(deltaNew, 0, sigma);
        path.PutInBox(deltaNew);
        beadB(i)->r = rBarNew + deltaNew;

        // Get sampling probs
        gaussProdOld = 1.;
        gaussProdNew = 1.;
        for (int iD=0; iD<path.nD; iD++) {
          gaussSumOld = 0.;
          gaussSumNew = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            distOld = deltaOld(iD) + (RealType)image*path.L;
            distNew = deltaNew(iD) + (RealType)image*path.L;
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
    RealType oldAction = 0.;
    RealType newAction = 0.;
    for (int iAction=0; iAction<actionList.size(); ++iAction) {
      // Old action
      path.SetMode(0);
      oldAction += actionList[iAction]->GetAction(bead0, bead1, particles, iLevel);

      // New action
      path.SetMode(1);
      newAction += actionList[iAction]->GetAction(bead0, bead1, particles, iLevel);
    }

    // Calculate acceptance ratio
    RealType logSampleRatio = -newLogSampleProb + oldLogSampleProb;
    RealType currActionChange = newAction - oldAction;
    RealType logAcceptProb = logSampleRatio - currActionChange + prevActionChange;

    // Metropolis step
    if (logAcceptProb < log(rng.unifRand()))
      return 0;

    prevActionChange = currActionChange;
  }

  // Have to check this if neighborhood has changed
  if (nPermPart < nPart) {
    // Construct Permutation Table
    path.SetMode(1);
    RealType permTot1 = constructPermTable();
  
    // Decide whether or not to accept whole bisection
    if ((oldCycleWeight*permTot0/permTot1) < rng.unifRand())
      return 0;
  }

  return 1;
}

// Construct permutation table and probabilities
RealType PermBisect::constructPermTable()
{
  // Update t
  updatePermTable();

  // Run through permatation types
  RealType totalWeight = 0.;
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
    b0(iP) = path.bead(iP+offset,bead0);
    b1(iP) = path.GetNextBead(b0(iP),nBisectBeads);
  }

  // Construct t table
  Tvector dr_ij(path.nD), dr_ii(path.nD);
  RealType exponent;
  for (unsigned int i=0; i<nPart; i++) {
    path.Dr(b0(i), b1(i), dr_ii);
    for (unsigned int j=0; j<nPart; j++) {
      path.Dr(b0(i), b1(j), dr_ij);
      exponent = (-dot(dr_ij, dr_ij) + dot(dr_ii, dr_ii))*i4LambdaTauNBisectBeads;
      t(i,j) = exp(exponent);
    }
  }

}

int PermBisect::selectCycle(RealType permTot)
{
  RealType x = rng.unifRand(0.,permTot);
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
//  for (unsigned int iP=offset; iP<offset+nPart; iP++) {
//    b = path.bead(iP,0);
//    for (unsigned int iB=0; iB<path.nBead; iB++) {
//      path.bead(iP,iB) = b;
//      path.bead(iP,iB)->p = iP;
//      b = b->next;
//    }
//  }
  for (unsigned int iP=offset; iP<offset+nPart; iP++) {
    b = path.bead(iP,path.beadLoop(bead1-1));
    for (unsigned int iB=path.beadLoop(bead1-1); iB<path.nBead; iB++) {
      path.bead(iP,iB) = b;
      path.bead(iP,iB)->p = iP;
      b = b->next;
    }
  }

  //for (unsigned int iP=offset; iP<offset+nPart; iP++) {
  //  for (unsigned int iB=0; iB<path.nBead; iB++) {
  //    cout << iP << " " << iB << "   " << path.bead(iP,iB)->prev->p << " " << path.bead(iP,iB)->p << " " << path.bead(iP,iB)->next->p << "   " << path.bead(iP,iB)->prev->b << " " << path.bead(iP,iB)->b << " " << path.bead(iP,iB)->next->b << endl;
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
