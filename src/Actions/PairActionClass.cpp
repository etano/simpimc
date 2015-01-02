#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder",-1);
  speciesA = in.getAttribute<string>("speciesA");
  speciesList.push_back(speciesA);
  speciesB = in.getAttribute<string>("speciesB");
  speciesList.push_back(speciesB);
  cout << "Setting up pair action between " << speciesA << " and " << speciesB << "..." << endl;
  maxLevel = in.getAttribute<int>("maxLevel",0);
  useLongRange = in.getAttribute<int>("useLongRange",0);
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
    vec<double> dr(path.nD);
    if (iSpeciesA == iSpeciesB) {
      for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
        for (int jP=iP+1; jP<path.speciesList[iSpeciesA]->nPart; ++jP) {
          for (int iB=0; iB<path.nBead; iB+=1) {
            int jB = iB + 1;
            path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesA,jP,iB),dr);
            double rMag = mag(dr);
            path.Dr(path(iSpeciesA,iP,jB),path(iSpeciesA,jP,jB),dr);
            double rPMag = mag(dr);
            tot += CalcV(rMag,rPMag,0);
          }
        }
      }
    } else {
      for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (int jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
          for (int iB=0; iB<path.nBead; iB+=1) {
            int jB = iB + 1;
            path.Dr(path(iSpeciesA,iP,iB),path(iSpeciesB,jP,iB),dr);
            double rMag = mag(dr);
            path.Dr(path(iSpeciesA,iP,jB),path(iSpeciesB,jP,jB),dr);
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
        for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
          #pragma omp for collapse(2) reduction(+:tot)
          for (int jP=iP+1; jP<path.speciesList[iSpeciesA]->nPart; ++jP) {
            for (int iB=0; iB<path.nBead; ++iB) {
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
      for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (int jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
          for (int iB=0; iB<path.nBead; ++iB) {
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

double PairAction::GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{

  if (level > maxLevel || isConstant || iSpeciesA < 0 || iSpeciesB < 0)
    return 0.;

  // Make sure particles are of species A or B and organize them accordingly
  vector<int> particlesA, particlesB;
  int nA(0), nB(0);
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpeciesA) {
      particlesA.push_back(iP);
      nA++;
    }
    else if (iS == iSpeciesB) {
      particlesB.push_back(iP);
      nB++;
    }
  }
  if (nA==0)
    if ((iSpeciesA==iSpeciesB) || (nB==0))
      return 0.;

  // Make vectors of other particles of species A and B
  vector<int> otherParticlesA;
  for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
    if (find(particlesA.begin(), particlesA.end(), iP)==particlesA.end())
      otherParticlesA.push_back(iP);
  }
  vector<int> otherParticlesB;
  for (int iP=0; iP<path.speciesList[iSpeciesB]->nPart; ++iP) {
    if (find(particlesB.begin(), particlesB.end(), iP)==particlesB.end())
      otherParticlesB.push_back(iP);
  }

  // Sum up contributing terms
  int skip = 1<<level;
  double levelTau = skip*path.tau;
  double tot = 0.;

  
  // Homologous
  if (iSpeciesA == iSpeciesB) {
    // Loop over beads
    for (int iB=b0; iB<b1; iB+=skip) {
      int jB = iB+skip;
      // Loop over A particles with other A particles
      for (auto& p: particlesA) {
        for (auto& q: otherParticlesA) {
          double rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesA,p,q,rMag,rPMag,rrPMag);
          double U = CalcU(rMag,rPMag,rrPMag,level);
          tot += U;
        }
      }
      // Loop over A particles with A particles
      for (int p=0; p<nA-1; ++p) {
        for (int q=p+1; q<nA; ++q) {
          double rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,particlesA[p],particlesA[q],rMag,rPMag,rrPMag);
          double U = CalcU(rMag,rPMag,rrPMag,level);
          tot += U;
        }
      }
    }
  // Heterologous
  } else {
    // Loop over beads
    for (int iB=b0; iB<b1; iB+=skip) {
      int jB = iB+skip;
      // Loop over A particles with other B particles
      for (auto& p: particlesA) {
        for (auto& q: otherParticlesB) {
          double rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,p,q,rMag,rPMag,rrPMag);
          double U = CalcU(rMag,rPMag,rrPMag,level);
          tot += U;
        }
      }
      // Loop other A particles with B particles
      for (auto& q: otherParticlesA) {
        for (auto& p: particlesB) {
          double rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,q,p,rMag,rPMag,rrPMag);
          double U = CalcU(rMag,rPMag,rrPMag,level);
          tot += U;
        }
      }
      // Loop over A particles with B particles
      for (auto& p: particlesA) {
        for (auto& q: particlesB) {
          double rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iSpeciesA,iSpeciesB,p,q,rMag,rPMag,rrPMag);
          double U = CalcU(rMag,rPMag,rrPMag,level);
          tot += U;
        }
      }
    }
  }

  // Add in long range part
  if (useLongRange) { // fixme: currently this assumes level = 0
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
