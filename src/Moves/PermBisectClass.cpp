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
  epsilon = in.getAttribute<RealType>("epsilon",1.e-100);
  logEpsilon = log(epsilon);

  // Set species things
  path.GetSpeciesInfo(species,iSpecies,offset);

  // Initialize constant cofactor
  nBisectBeads = 1<<nLevel; // Number of beads in bisection
  lambda = path.speciesList[iSpecies]->lambda;
  i4LambdaTauNBisectBeads = 1./(4.*lambda*path.tau*nBisectBeads);

  // Initiate permutation table
  nPart = path.speciesList[iSpecies]->nPart;
  t.zeros(nPart,nPart);
  nPermPart = 3; // fixme: make more general
  int Nchoose3 = nPart*(nPart-1)*(nPart-2)/6;
  all_cycles.set_size(path.nPermType * Nchoose3);
  BuildCycles();

  // Initiate acceptance ratio counters
  permAttempt.zeros(path.nPermType);
  permAccept.zeros(path.nPermType);
}

void PermBisect::MakeMove()
{
  nAccept += DoPermBisect();
  nAttempt++;
}

int PermBisect::DoPermBisect()
{
  unsigned int bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  unsigned int bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  bool rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle

  // Set up permutation
  cycles.clear();
  RealType permTot0 = constructPermTable(bead0,nBisectBeads); // Permutation weight table
  int cycleIndex = selectCycle(permTot0);
  Cycle* c = cycles[cycleIndex];

  // Set up pointers
  vector<int> particles;
  field<Bead*> beadI(nPermPart), beadFm1(nPermPart), beadF(nPermPart);
  for (unsigned int i=0; i<nPermPart; i++) {
    beadI(i) = path.bead(c->part(i),bead0);
    beadFm1(i) = beadI(i)->nextB(nBisectBeads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(c->part(i));
  }

  // Permute particles
  permuteBeads(beadFm1, beadF, c);

  // Set final bisection beads and store
  field<Bead*> beadA(nPermPart);
  affBeads.clear();
  for (unsigned int i=0; i<nPermPart; i++) {
    for(beadA(i) = beadI(i); beadA(i) != beadF(i); beadA(i) = beadA(i) -> next)
      affBeads.push_back(beadA(i));
  }

  // Perform the bisection (move exactly through kinetic action)
  field<Bead*> beadB(nPermPart), beadC(nPermPart);
  RealType prevActionChange = 0.;
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
    for (unsigned int i=0; i<nPermPart; i++) {
      beadA(i) = beadI(i);
      while(beadA(i) != beadF(i)) {
        // Set beads
        beadB(i) = beadA(i)->nextB(skip);
        beadC(i) = beadB(i)->nextB(skip);

        // Old sampling
        path.SetMode(0);
        path.RBar(beadC(i), beadA(i), rBarOld);
        path.Dr(beadB(i), rBarOld, deltaOld);

        // New sampling
        path.SetMode(1);
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
            gaussSumOld += exp(-0.5*distOld*distOld/sigma2);
            gaussSumNew += exp(-0.5*distNew*distNew/sigma2);
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
    for (unsigned int iB=bead0; iB<bead1; iB+=2*skip) {
      for (int iAction=0; iAction<actionList.size(); ++iAction) {
        // Old action
        path.SetMode(0);
        oldAction += actionList[iAction]->GetAction(iB, iB+2*skip, particles, iLevel);

        // New action
        path.SetMode(1);
        newAction += actionList[iAction]->GetAction(iB, iB+2*skip, particles, iLevel);
      }
    }

    // Calculate acceptance ratio
    RealType logSampleRatio = -newLogSampleProb + oldLogSampleProb;
    RealType currActionChange = newAction - oldAction;
    RealType logAcceptProb = logSampleRatio - currActionChange + prevActionChange;

    // Metropolis step
    if (logAcceptProb < log(rng.unifRand())) { // Reject if true
      // Restore things
      for (unsigned int iP=0; iP<particles.size(); iP++) {
        path.bead(particles[iP],path.beadLoop(bead1)) -> restorePrev();
        path.bead(particles[iP],path.beadLoop(bead1-1)) -> restoreNext();
      }
      assignParticleLabels();
      path.restoreR(affBeads);
      path.restoreRhoK(affBeads);

      return 0;
    }

    prevActionChange = currActionChange;
  }

  // fixme: Decide whether or not it is necessary to have his term
  //        It seems it doesn't matter (with all gammas=1)
  //// Construct Permutation Table
  //path.SetMode(1);
  //RealType permTot1 = constructPermTable(bead0,nBisectBeads);

  //// Decide whether or not to accept whole bisection
  //if ((permTot0/permTot1) < rng.unifRand())  {
  //  // Restore things
  //  for (unsigned int iP=0; iP<particles.size(); iP++) {
  //    path.bead(particles[iP],path.beadLoop(bead1)) -> restorePrev();
  //    path.bead(particles[iP],path.beadLoop(bead1-1)) -> restoreNext();
  //  }
  //  assignParticleLabels();
  //  path.restoreR(affBeads);
  //  path.restoreRhoK(affBeads);

  //  return 0;
  //}

  // Accept move, so store things
  for (unsigned int iP=0; iP<particles.size(); iP++) {
    path.bead(particles[iP],path.beadLoop(bead1)) -> storePrev();
    path.bead(particles[iP],path.beadLoop(bead1-1)) -> storeNext();
  }
  assignParticleLabels();
  path.storeR(affBeads);
  path.storeRhoK(affBeads);

  // Increment permutation counter
  permAccept(c->type) += 1;

  return 1;
}

// Construct all cycles of 3 Particle Exchanges
void PermBisect::BuildCycles()
{
  // Run through permatation types
  int permIndex = 0;
  for (unsigned int permType=0; permType<path.nPermType; permType++) {
    for (unsigned int i=0; i<nPart-2; i++) {
      for (unsigned int j=i+1; j<nPart-1; j++) {
        for (unsigned int k=j+1; k<nPart; k++) {
          Cycle& c = all_cycles(permIndex);
          c.type = permType;
          c.index = permIndex;
          c.part = {i, j, k};
          setPerm(c);
          permIndex += 1;
        }
      }
    }
  }

}

// All 3 particle exchanges (first 1-2 for fermions and bosons, last 3-5 for bosons only, default is identity)
void PermBisect::setPerm(Cycle& c)
{
  // Permutations
  switch (c.type) {
    case 1:
      c.perm = {2, 0, 1};
      c.iPerm = {1, 2, 0};
      break;
    case 2:
      c.perm = {1, 2, 0};
      c.iPerm = {2, 0, 1};
      break;
    case 3:
      c.perm = {1, 0, 2};
      c.iPerm = {1, 0, 2};
      break;
    case 4:
      c.perm = {0, 2, 1};
      c.iPerm = {0, 2, 1};
      break;
    case 5:
      c.perm = {2, 1, 0};
      c.iPerm = {2, 1, 0};
      break;
    default:
      c.perm = {0, 1, 2};
      c.iPerm = {0, 1, 2};
      break;
  }
}

// Construct permutation table and probabilities
RealType PermBisect::constructPermTable(const int bead0, const int nBisectBeads)
{
  // Set initial and final beads
  field<Bead*> b0(nPart), b1(nPart);
  for (unsigned int iP=0; iP<nPart; iP++) {
    b0(iP) = path.bead(iP+offset,bead0);
    b1(iP) = b0(iP) -> nextB(nBisectBeads);
  }

  // Construct t table
  Tvector dr(path.nD);
  RealType exponent;
  for (unsigned int i=0; i<nPart; i++) {
    for (unsigned int j=0; j<nPart; j++) {
      path.Dr(b0(i), b1(j), dr);
      exponent = -dot(dr, dr)*i4LambdaTauNBisectBeads;
      if (exponent > logEpsilon)
        t(i,j) = exp(exponent);
      else
        t(i,j) = 0.;
    }
  }

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
  permAttempt(c->type) += 1;

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
  for (unsigned int iP=offset; iP<offset+nPart; iP++) {
    b = path.bead(iP,0);
    for (unsigned int iB=0; iB<path.nBead; iB++) {
      path.bead(iP,iB) = b;
      path.bead(iP,iB)->p = iP;
      b = b->next;
    }
  }
}

void PermBisect::Write()
{
  // Write
  if (firstTime) {
    out.Write("/Moves/"+name+"/nPermType", path.nPermType);
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
