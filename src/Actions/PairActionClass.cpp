#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<uint>("nOrder",0);
  speciesA = in.getAttribute<string>("speciesA");
  speciesList.push_back(speciesA);
  speciesB = in.getAttribute<string>("speciesB");
  speciesList.push_back(speciesB);
  cout << "Setting up pair action between " << speciesA << " and " << speciesB << "..." << endl;
  maxLevel = in.getAttribute<uint>("maxLevel",0);
  useLongRange = in.getAttribute<bool>("useLongRange",0);
  if (useLongRange) {
    kCut = in.getAttribute<double>("kCut",path.kC);
    path.SetupKs(kCut);
  }
  path.GetSpeciesInfo(speciesA,iSpeciesA);
  path.GetSpeciesInfo(speciesB,iSpeciesB);
  if (iSpeciesA >= 0 && iSpeciesB >= 0) {
    isConstant = ((iSpeciesA == iSpeciesB)
                 && (path.speciesList[iSpeciesA]->nPart == 1
                    || path.speciesList[iSpeciesA]->lambda == 0.));
    isFirstTime = true;

    string fileName = in.getAttribute<string>("file");
    ReadFile(fileName);

    // Write things to file
    out.Write("Actions/"+name+"/file", fileName);
    out.Write("Actions/"+name+"/nImages", nImages);
    out.Write("Actions/"+name+"/nOrder", nOrder);
    out.Write("Actions/"+name+"/speciesA", speciesA);
    out.Write("Actions/"+name+"/speciesB", speciesB);
    out.Write("Actions/"+name+"/maxLevel", maxLevel);
    out.Write("Actions/"+name+"/useLongRange", useLongRange);
    if (useLongRange)
      out.Write("Actions/"+name+"/kCut", kCut);
  } else {
    isConstant = true;
    isFirstTime = false;
    dUdBConstant = 0.;
    VConstant = 0.;
  }

}

double PairAction::Potential()
{
  if (isConstant && !isFirstTime)
    return VConstant;
  else {
    double tot = 0.;
    if (iSpeciesA == iSpeciesB) {
      for (uint iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
        for (uint jP=iP+1; jP<path.speciesList[iSpeciesA]->nPart; ++jP) {
          for (uint iB=0; iB<path.nBead; iB+=1) {
            uint jB = iB + 1;
            vec<double> dr(path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesA,jP,iB)));
            double rMag = mag(dr);
            dr = path.Dr(path(iSpeciesA,iP,jB),path(iSpeciesA,jP,jB));
            double rPMag = mag(dr);
            tot += CalcV(rMag,rPMag,0);
          }
        }
      }
    } else {
      for (uint iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (uint jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
          for (uint iB=0; iB<path.nBead; iB+=1) {
            uint jB = iB + 1;
            vec<double> dr(path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesB,jP,iB)));
            double rMag = mag(dr);
            dr = path.Dr(path(iSpeciesA,iP,jB),path(iSpeciesB,jP,jB));
            double rPMag = mag(dr);
            tot += CalcV(rMag,rPMag,0);
          }
        }
      }
    }

    if (useLongRange)
      tot += CalcVLong();

    if (isFirstTime) {
      isFirstTime = false;
      VConstant = tot;
    }

    return tot;
  }
}

double PairAction::DActionDBeta()
{
  if (isConstant && !isFirstTime)
    return dUdBConstant;
  else {
    double tot = 0.;
    if (iSpeciesA == iSpeciesB) {
      #pragma omp parallel
      {
        for (uint iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
          #pragma omp for collapse(2) reduction(+:tot)
          for (uint jP=iP+1; jP<path.speciesList[iSpeciesA]->nPart; ++jP) {
            for (uint iB=0; iB<path.nBead; ++iB) {
              double rMag, rPMag, rrPMag;
              path.DrDrPDrrP(iB,iB+1,iSpeciesA,iSpeciesA,iP,jP,rMag,rPMag,rrPMag);
              double dUdB = CalcdUdBeta(rMag,rPMag,rrPMag,0);
              tot += dUdB;
            }
          }
        }
      }
    } else {
      #pragma omp parallel for collapse(3) reduction(+:tot)
      for (uint iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (uint jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
          for (uint iB=0; iB<path.nBead; ++iB) {
            double rMag, rPMag, rrPMag;
            path.DrDrPDrrP(iB,iB+1,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
            double dUdB = CalcdUdBeta(rMag,rPMag,rrPMag,0);
            tot += dUdB;
          }
        }
      }
    }

    if (useLongRange)
      tot += CalcdUdBetaLong();

    if (isFirstTime) {
      isFirstTime = false;
      dUdBConstant = tot;
    }

    return tot;
  }
}

void PairAction::GenerateParticlePairs(const vector<pair<uint,uint> > &particles, vector<uint> &particlesA, vector<uint> &particlesB, vector< pair<uint,uint> > &particlePairs)
{
  // Make sure particles are of species A or B and organize them accordingly
  uint nA(0), nB(0);
  for (auto& p: particles) {
    uint iS = p.first;
    uint iP = p.second;
    if (iS == iSpeciesA) {
      particlesA.push_back(iP);
      nA++;
    } else if (iS == iSpeciesB) {
      particlesB.push_back(iP);
      nB++;
    }
  }
  if (nA==0)
    if ((iSpeciesA==iSpeciesB) || (nB==0))
      return;

  // Make vectors of other particles of species A
  vector<uint> otherParticlesA;
  for (uint iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
    if (find(particlesA.begin(), particlesA.end(), iP)==particlesA.end())
      otherParticlesA.push_back(iP);
  }
  // Make vectors of other particles of species B
  vector<uint> otherParticlesB;
  for (uint iP=0; iP<path.speciesList[iSpeciesB]->nPart; ++iP) {
    if (find(particlesB.begin(), particlesB.end(), iP)==particlesB.end())
      otherParticlesB.push_back(iP);
  }

  // Homologous
  if (iSpeciesA == iSpeciesB) {
    // Loop over A particles with other A particles
    for (auto& p: particlesA)
      for (auto& q: otherParticlesA)
        particlePairs.push_back(std::make_pair(p,q));
    // Loop over A particles with A particles
    for (uint p=0; p<nA-1; ++p)
      for (uint q=p+1; q<nA; ++q)
        particlePairs.push_back(std::make_pair(particlesA[p],particlesA[q]));
  // Heterologous
  } else {
    // Loop over A particles with other B particles
    for (auto& p: particlesA)
      for (auto& q: otherParticlesB)
        particlePairs.push_back(std::make_pair(p,q));
    // Loop other A particles with B particles
    for (auto& p: otherParticlesA)
      for (auto& q: particlesB)
        particlePairs.push_back(std::make_pair(p,q));
    // Loop over A particles with B particles
    for (auto& p: particlesA)
      for (auto& q: particlesB)
        particlePairs.push_back(std::make_pair(p,q));
  }

}

double PairAction::GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> > &particles, const uint level)
{
  // Return zero if not relevant
  if (level > maxLevel || isConstant || iSpeciesA < 0 || iSpeciesB < 0)
    return 0.;

  // Generate particle pairs
  vector<uint> particlesA, particlesB;
  vector< pair<uint,uint> > particlePairs;
  GenerateParticlePairs(particles, particlesA, particlesB, particlePairs);
  if (particlePairs.size() == 0)
    return 0.;

  // Sum up contributing terms
  uint skip = 1<<level;
  double tot = 0.;
  #pragma omp parallel for reduction(+:tot)
  for (uint iB=b0; iB<b1; iB+=skip) {
    uint jB = iB+skip;
    uint kB = iB-skip;
    if (b0 == 0)
      kB = path.nBead-1;
    for (auto& particlePair: particlePairs) {
      double rMag, rPMag, rrPMag;
      path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,particlePair.first,particlePair.second,rMag,rPMag,rrPMag);
      tot += CalcU(rMag,rPMag,rrPMag,level);
    }
  }

  // Add in long range part
  if (useLongRange) { // FIXME: currently this assumes level = 0
    if (path.speciesList[iSpeciesA]->needUpdateRhoK && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, iSpeciesA, particlesA, level);
      path.speciesList[iSpeciesA]->needUpdateRhoK = false;
    }
    if (path.speciesList[iSpeciesB]->needUpdateRhoK && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, iSpeciesB, particlesB, level);
      path.speciesList[iSpeciesB]->needUpdateRhoK = false;
    }
    tot += CalcULong(b0, b1, level);
  }

  return tot;
}

vec<double> PairAction::GetActionGradient(const uint b0, const uint b1, const vector<pair<uint,uint> > &particles, const uint level)
{
  // Return zero if not relevant
  vec<double> zero_vec;
  zero_vec.zeros(path.nD);
  if (level > maxLevel || isConstant || iSpeciesA < 0 || iSpeciesB < 0)
    return zero_vec;

  // Generate particle pairs
  vector<uint> particlesA, particlesB;
  vector< pair<uint,uint> > particlePairs;
  GenerateParticlePairs(particles, particlesA, particlesB, particlePairs);
  if (particlePairs.size() == 0)
    return zero_vec;

  // Sum up contributing terms
  uint skip = 1<<level;
  vec<double> tot(zero_vec);
  for (uint iB=b0; iB<b1; iB+=skip) {
    uint jB = iB+skip;
    uint kB = iB-skip;
    if (iB == 0) // FIXME: This doesn't depend on the level
      kB = path.nBead-1;
    for (auto& particlePair: particlePairs)
      tot += CalcGradientU(iB,jB,particlePair.first,particlePair.second,level)
             + CalcGradientU(iB,kB,particlePair.first,particlePair.second,level);
  }

  // FIXME: Ignoring long range part for now
  //// Add in long range part
  //if (useLongRange) { // fixme: currently this assumes level = 0
  //  if (path.speciesList[iSpeciesA]->needUpdateRhoK && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, iSpeciesA, particlesA, level);
  //    path.speciesList[iSpeciesA]->needUpdateRhoK = false;
  //  }
  //  if (path.speciesList[iSpeciesB]->needUpdateRhoK && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, iSpeciesB, particlesB, level);
  //    path.speciesList[iSpeciesB]->needUpdateRhoK = false;
  //  }
  //  tot += CalcGradientULong(b0, b1, level);
  //}

  return tot;
}

double PairAction::GetActionLaplacian(const uint b0, const uint b1, const vector<pair<uint,uint> > &particles, const uint level)
{
  // Return zero if not relevant
  if (level > maxLevel || isConstant || iSpeciesA < 0 || iSpeciesB < 0)
    return 0.;

  // Generate particle pairs
  vector<uint> particlesA, particlesB;
  vector< pair<uint,uint> > particlePairs;
  GenerateParticlePairs(particles, particlesA, particlesB, particlePairs);
  if (particlePairs.size() == 0)
    return 0.;

  // Sum up contributing terms
  uint skip = 1<<level;
  double tot = 0.;
  for (uint iB=b0; iB<b1; iB+=skip) {
    uint jB = iB+skip;
    uint kB = iB-skip;
    if (iB == 0) // FIXME: This doesn't depend on the level
      kB = path.nBead-1;
    for (auto& particlePair: particlePairs)
      tot += CalcLaplacianU(iB,jB,particlePair.first,particlePair.second,level)
             + CalcLaplacianU(iB,kB,particlePair.first,particlePair.second,level);
  }

  // FIXME: Ignoring long range part for now
  //// Add in long range part
  //if (useLongRange) { // fixme: currently this assumes level = 0
  //  if (path.speciesList[iSpeciesA]->needUpdateRhoK && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, iSpeciesA, particlesA, level);
  //    path.speciesList[iSpeciesA]->needUpdateRhoK = false;
  //  }
  //  if (path.speciesList[iSpeciesB]->needUpdateRhoK && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, iSpeciesB, particlesB, level);
  //    path.speciesList[iSpeciesB]->needUpdateRhoK = false;
  //  }
  //  tot += CalcLaplacianULong(b0, b1, level);
  //}

  return tot;
}

vec<double> PairAction::CalcGradientU(const uint iB, const uint jB, const uint iP, const uint jP, const uint level)
{
  // Store original position for particle i
  std::shared_ptr<Bead> b(path(iSpeciesA,iP,iB));
  vec<double> r0(b->r);

  // Numerical tolerance
  double eps = 1.e-7;

  // Calculate numerical gradient
  double rMag, rPMag, rrPMag;
  vec<double> tot(path.nD);
  for (uint iD=0; iD<path.nD; ++iD) {
    b->r(iD) = r0(iD) + eps;
    path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
    double f1 = CalcU(rMag,rPMag,rrPMag,level);
    b->r(iD) = r0(iD) - eps;
    path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
    double f2 = CalcU(rMag,rPMag,rrPMag,level);
    tot(iD) = (f1-f2)/(2.*eps);
    b->r(iD) = r0(iD);
  }

  return tot;
}

double PairAction::CalcLaplacianU(const uint iB, const uint jB, const uint iP, const uint jP, const uint level)
{
  // Store original position for particle i
  std::shared_ptr<Bead> b(path(iSpeciesA,iP,iB));
  vec<double> r0(b->r);

  // Numerical tolerance
  double eps = 1.e-7;

  // Calculate numerical gradient
  double rMag, rPMag, rrPMag;
  double tot = 0.;
  for (uint iD=0; iD<path.nD; ++iD) {
    path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
    double f0 = CalcU(rMag,rPMag,rrPMag,level);
    b->r(iD) = r0(iD) + eps;
    path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
    double f1 = CalcU(rMag,rPMag,rrPMag,level);
    b->r(iD) = r0(iD) - eps;
    path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,iP,jP,rMag,rPMag,rrPMag);
    double f2 = CalcU(rMag,rPMag,rrPMag,level);
    tot += (f1+f2-2.*f0)/(eps*eps);
    b->r(iD) = r0(iD);
  }

  return tot;
}

double PairAction::ImportanceWeight()
{
  return isImportanceWeight ? exp(DActionDBeta()/path.nBead) : 1.;
}

void PairAction::Write() {}

void PairAction::Accept() // fixme: will accept even when unaffected
{
  if (useLongRange && iSpeciesA >= 0 && iSpeciesB >= 0) {
    path.speciesList[iSpeciesA]->needUpdateRhoK = true;
    path.speciesList[iSpeciesB]->needUpdateRhoK = true;
  }

}

void PairAction::Reject() // fixme: will accept even when unaffected, why is this true?
{
  if (useLongRange && iSpeciesA >= 0 && iSpeciesB >= 0) {
    path.speciesList[iSpeciesA]->needUpdateRhoK = true;
    path.speciesList[iSpeciesB]->needUpdateRhoK = true;
  }

}
