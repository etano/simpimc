#include "PermBisectIterativeClass.h"

void PermBisectIterative::Init(Input &in)
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
  nPart = path.speciesList[iSpecies]->nPart;

  // Initialize constant cofactor
  nBisectBeads = 1<<nLevel; // Number of beads in bisection
  lambda = path.speciesList[iSpecies]->lambda;
  i4LambdaTauNBisectBeads = 1./(4.*lambda*path.tau*nBisectBeads);

  // Initiate permutation table
  t.zeros(nPart,nPart);

  // Initiate acceptance ratio counters
  permAttempt.zeros(nPart);
  permAccept.zeros(nPart);
}

void PermBisectIterative::Accept()
{
  permAttempt(nPermPart-1) += 1;
  permAccept(nPermPart-1) += 1;

  // Change sign weight for fermions
  if (!(nPermPart%2) && path.speciesList[iSpecies]->fermi)
    path.sign *= -1;

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

  // Call reject for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Accept();
}

void PermBisectIterative::Reject()
{
  // No need to do some things if bisection isn't attempted
  if (nPermPart > 0) {
    permAttempt(nPermPart-1) += 1;

    // Restore things
    for (unsigned int iP=offset; iP<offset+nPart; iP++) { // todo: can make this more efficient by only restoring touched particles
      path.bead(iP,path.beadLoop(bead1)) -> restorePrev();
      path.bead(iP,path.beadLoop(bead1-1)) -> restoreNext();
    }
    path.restoreR(affBeads);
    path.restoreRhoKP(affBeads);
    for (int iB=bead0; iB<=bead1; ++iB)
      path.restoreRhoK(iB,iSpecies);
  }

  // Call reject for each action
  for (int iAction=0; iAction<actionList.size(); ++iAction)
    actionList[iAction]->Reject();
}

// Perform the permuting bisection
int PermBisectIterative::Attempt()
{
  bead0 = rng.unifRand(path.nBead) - 1;  // Pick first bead at random
  bead1 = bead0 + nBisectBeads; // Set last bead in bisection
  rollOver = bead1 > (path.nBead-1);  // See if bisection overflows to next particle

  // Set up permutation
  Cycle c;
  nPermPart = 0; // reset to indicate if bisection is atempted or not
  if (!selectCycleIterative(c))
    return 0; // do not attempt bisection since permutation not accepted
  nPermPart = c.part.size();

  // Set up pointers
  vector<int> particles;
  field<Bead*> beadI(nPermPart), beadFm1(nPermPart), beadF(nPermPart);
  for (unsigned int i=0; i<nPermPart; i++) {
    beadI(i) = path.bead(c.part(i),bead0);
    beadFm1(i) = beadI(i)->nextB(nBisectBeads-1);
    beadF(i) = beadFm1(i)->next;
    particles.push_back(c.part(i));
  }
  if (rollOver) {
    for (int i=0; i<nPermPart; ++i)
      particles.push_back(beadF(i)->p);
    sort(particles.begin(), particles.end());
    particles.erase(unique(particles.begin(), particles.end()), particles.end());
  }

  // Permute particles
  permuteBeads(beadFm1, beadF, c);

  // Note affected beads
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

  return 1;
}

void PermBisectIterative::updatePermTable()
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
  Tmatrix t_c = t;

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
    t_c(p,p0) = t(p,p0);

    // Calculate row total
    RealType Q_p = 0.;
    RealType Q_p_c = 0.;
    for (int i=0; i<nPart; ++i) {
      Q_p += t(p,i);
      Q_p_c += t_c(p,i);
    }

    // Decide whether or not to continue
    if ((Q_p_c/Q_p) < rng.unifRand())
      return 0;

    // Select next particle with bisective search
    RealType x = rng.unifRand(0.,Q_p_c);
    RealType t_Q = 0.;
    for (int i=0; i<nPart; ++i) { // fixme: not doing bisection
      t_Q += t_c(p,i);
      if (t_Q > x) {
        p = i;
        break;
      }
    }

  } while (p != p0);

  // Disallow even permutations from closing for fixed-node calculations
  if (path.speciesList[iSpecies]->fermi && path.speciesList[iSpecies]->fixedNode && !(ps.size()%2)) {
    nPermPart = ps.size();
    return 0;
  }

  // Set particles
  int nPerm = ps.size();
  c.part = ps;
  for (int i=0; i<nPerm; ++i)
    c.part(i) += offset;

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
void PermBisectIterative::permuteBeads(field<Bead*>& b0, field<Bead*>& b1, Cycle& c)
{
  // Execute the permutation
  int nPerm = c.part.size();
  for (unsigned int i=0; i<nPerm; i++) {
    b0(i)->next = b1(c.perm(i));
  }
  for (unsigned int i=0; i<nPerm; i++)
    b1(i)->prev = b0(c.iPerm(i));
  for (unsigned int i=0; i<nPerm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisectIterative::assignParticleLabels()
{
  Bead *b;
  if (path.nBead-path.beadLoop(bead1-1) < path.beadLoop(bead1+1)) {
    for (unsigned int iP=offset; iP<offset+nPart; iP++) {
      b = path.bead(iP,path.beadLoop(bead1-1));
      for (unsigned int iB=path.beadLoop(bead1-1); iB<path.nBead; iB++) {
        path.bead(iP,iB) = b;
        path.bead(iP,iB)->p = iP;
        b = b->next;
      }
    }
  } else {
    for (unsigned int iP=offset; iP<offset+nPart; iP++) {
      b = path.bead(iP,path.beadLoop(bead1+1));
      for (unsigned int iB=path.beadLoop(bead1+1); iB>0; iB--) {
        path.bead(iP,iB) = b;
        path.bead(iP,iB)->p = iP;
        b = b->prev;
      }
    }
  }
}

void PermBisectIterative::Write()
{
  // Write
  if (firstTime) {
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
