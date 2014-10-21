#include "PairActionClass.h"

void PairAction::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  nOrder = in.getAttribute<int>("nOrder");
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  cout << "Setting up pair action between " << speciesA << " and " << speciesB << "..." << endl;
  maxLevel = in.getAttribute<int>("maxLevel",0);
  useLongRange = in.getAttribute<int>("useLongRange",0);
  if (useLongRange) {
    kCut = in.getAttribute<RealType>("kCut");
    path.SetupKs(kCut);
  }
  GetOffset(speciesA,iSpeciesA,offsetA);
  GetOffset(speciesB,iSpeciesB,offsetB);
  cout << speciesA << " " << speciesB << " " << iSpeciesA << " " << iSpeciesB << endl;
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
    dUdB = 0.;
  }

}

RealType PairAction::DActionDBeta()
{
  if (isConstant && !isFirstTime)
    return dUdB;
  else {
    RealType tot = 0.;
    Tvector r(path.nD), rP(path.nD), rrP(path.nD);
    if (iSpeciesA == iSpeciesB) {
      for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart-1; ++iP) {
        for (int jP=iP+1; jP<offsetA+path.speciesList[iSpeciesA]->nPart; ++jP) {
          for (int iB=0; iB<path.nBead; iB+=1) {
            int jB = iB + 1;
            RealType rMag, rPMag, rrPMag;
            path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
            tot += CalcdUdBeta(rMag,rPMag,rrPMag,0);
          }
        }
      }
    } else {
      for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
        for (int jP=offsetB; jP<offsetB+path.speciesList[iSpeciesB]->nPart; ++jP) {
          for (int iB=0; iB<path.nBead; iB+=1) {
            int jB = iB + 1;
            RealType rMag, rPMag, rrPMag;
            path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
            tot += CalcdUdBeta(rMag,rPMag,rrPMag,0);
          }
        }
      }
    }

    if (useLongRange)
      tot += CalcdUdBetaLong();

    isFirstTime = false;
    dUdB = tot;

    return tot;
  }
}

RealType PairAction::GetAction(int b0, int b1, vector<int> &particles, int level)
{

  if (level > maxLevel || isConstant || iSpeciesA < 0 || iSpeciesB < 0)
    return 0.;

  // Make sure particles are of species A or B and organize them accordingly
  vector<int> particlesA, particlesB;
  int nA(0), nB(0);
  for (int p=0; p<particles.size(); ++p) {
    int iP = particles[p];
    int iS = path(iP,b0)->s;
    if (iS == iSpeciesA) {
      particlesA.push_back(iP);
      nA++;
    }
    else if (iS == iSpeciesB) {
      particlesB.push_back(iP);
      nB++;
    }
  }
  if (nA==0 && nB==0)
    return 0.;
  if (iSpeciesA == iSpeciesB && (nA==0))
    return 0.;

  // Make vectors of other particles of species A and B
  vector<int> otherParticlesA;
  for (int iP=offsetA; iP<offsetA+path.speciesList[iSpeciesA]->nPart; ++iP) {
    if (find(particlesA.begin(), particlesA.end(), iP)==particlesA.end())
      otherParticlesA.push_back(iP);
  }
  vector<int> otherParticlesB;
  for (int iP=offsetB; iP<offsetB+path.speciesList[iSpeciesB]->nPart; ++iP) {
    if (find(particlesB.begin(), particlesB.end(), iP)==particlesB.end())
      otherParticlesB.push_back(iP);
  }

  // Sum up contributing terms
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  Tvector r(path.nD), rP(path.nD), rrP(path.nD);

  // Homologous
  if (iSpeciesA == iSpeciesB) {
    // Loop over A particles with other A particles
    for (int p=0; p<nA; ++p) {
      int iP = particlesA[p];
      for (int q=0; q<otherParticlesA.size(); ++q) {
        int jP = otherParticlesA[q];
        for (int iB=b0; iB<b1; iB+=skip) {
          int jB = iB + skip;
          RealType rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
          tot += CalcU(rMag,rPMag,rrPMag,level);
        }
      }
    }
    // Loop over A particles with A particles
    for (int p=0; p<nA-1; ++p) {
      int iP = particlesA[p];
      for (int q=p+1; q<nA; ++q) {
        int jP = particlesA[q];
        for (int iB=b0; iB<b1; iB+=skip) {
          int jB = iB + skip;
          RealType rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
          tot += CalcU(rMag,rPMag,rrPMag,level);
        }
      }
    }
  // Heterologous
  } else {
    // Loop over A particles with other B particles
    for (int p=0; p<nA; ++p) {
      int iP = particlesA[p];
      for (int q=0; q<otherParticlesB.size(); ++q) {
        int jP = otherParticlesB[q];
        for (int iB=b0; iB<b1; iB+=skip) {
          int jB = iB + skip;
          RealType rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
          tot += CalcU(rMag,rPMag,rrPMag,level);
        }
      }
    }
    // Loop over B particles with other A particles
    for (int p=0; p<nB; ++p) {
      int iP = particlesB[p];
      for (int q=0; q<otherParticlesA.size(); ++q) {
        int jP = otherParticlesA[q];
        for (int iB=b0; iB<b1; iB+=skip) {
          int jB = iB + skip;
          RealType rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
          tot += CalcU(rMag,rPMag,rrPMag,level);
        }
      }
    }
    // Loop over A particles with B particles
    for (int p=0; p<nA; ++p) {
      int iP = particlesA[p];
      for (int q=0; q<nB; ++q) {
        int jP = particlesB[q];
        for (int iB=b0; iB<b1; iB+=skip) {
          int jB = iB + skip;
          RealType rMag, rPMag, rrPMag;
          path.DrDrPDrrP(iB,jB,iP,jP,rMag,rPMag,rrPMag,r,rP,rrP);
          tot += CalcU(rMag,rPMag,rrPMag,level);
        }
      }
    }
  }

  // Add in long range part
  if (useLongRange) { // fixme: currently this assumes level = 0
    if (path.speciesList[iSpeciesA]->needUpdateRhoK && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, particlesA, level);
      path.speciesList[iSpeciesA]->needUpdateRhoK = false;
    }
    if (path.speciesList[iSpeciesB]->needUpdateRhoK && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, particlesB, level);
      path.speciesList[iSpeciesB]->needUpdateRhoK = false;
    }
    tot += CalcULong(b0, b1, particles, level);
  }

  return tot;
}

void PairAction::Write() {}

void PairAction::Accept()
{
  if (useLongRange && iSpeciesA >= 0 && iSpeciesB >= 0) {
    path.speciesList[iSpeciesA]->needUpdateRhoK = true;
    path.speciesList[iSpeciesB]->needUpdateRhoK = true;
  }

}

void PairAction::Reject()
{
  if (useLongRange && iSpeciesA >= 0 && iSpeciesB >= 0) {
    path.speciesList[iSpeciesA]->needUpdateRhoK = true;
    path.speciesList[iSpeciesB]->needUpdateRhoK = true;
  }

}
