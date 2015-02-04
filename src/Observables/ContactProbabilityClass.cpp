#include "ContactProbabilityClass.h"

void ContactProbability::Init(Input &in)
{
  // Read in species info
  speciesA = in.getAttribute<string>("speciesA");
  speciesB = in.getAttribute<string>("speciesB");
  path.GetSpeciesInfo(speciesA, iSpeciesA);
  path.GetSpeciesInfo(speciesB, iSpeciesB);

  // Write things to file
  out.Write(prefix+"/speciesA", speciesA);
  out.Write(prefix+"/speciesB", speciesB);

  // Read in ZA
  ZA = in.getAttribute<int>("ZA");

  // Generate action list
  std::vector<std::string> speciesList;
  speciesList.push_back(speciesA);
  speciesList.push_back(speciesB);
  for (auto& action: fullActionList) {
    for (auto& sA: speciesList) {
      if (std::find(action->speciesList.begin(), action->speciesList.end(), sA) != action->speciesList.end()) {
        actionList.push_back(action);
        break;
      }
    }
  }

  Reset();
}

void ContactProbability::Reset()
{
  total = 0;
  nMeasure = 0;
}

// Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
void ContactProbability::Accumulate()
{
  path.SetMode(1);

  // Form particle pairs
  std::vector<std::vector<pair<int,int> > > particlePairs;
  if (iSpeciesA == iSpeciesB) { // Homogeneous
    for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart-1; ++iP) {
      for (int jP=iP+1; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
        std::vector<pair<int,int> > particles;
        particles.push_back(std::make_pair(iSpeciesA,iP));
        particles.push_back(std::make_pair(iSpeciesB,jP));
        particlePairs.push_back(particles);
      }
    }
  } else { // Homologous
    for (int iP=0; iP<path.speciesList[iSpeciesA]->nPart; ++iP) {
      for (int jP=0; jP<path.speciesList[iSpeciesB]->nPart; ++jP) {
        std::vector<pair<int,int> > particles;
        particles.push_back(std::make_pair(iSpeciesA,iP));
        particles.push_back(std::make_pair(iSpeciesB,jP));
        particlePairs.push_back(particles);
      }
    }
  }

  // Add up contact probability
  // FIXME: Currently only looking at origin
  for (auto& particles: particlePairs) {
    for (int iB=0; iB<path.nBead; ++iB) {
      // Set r's
      vec<double> RA = path(particles[0].first,particles[0].second,iB)->r;
      vec<double> ri = path(particles[1].first,particles[1].second,iB)->r;

      // Get differences
      vec<double> ri_RA(path.nD);
      path.Dr(ri, RA, ri_RA);
      double mag_ri_RA = mag(ri_RA);

      // Compute functions
      double g = 0.; // FIXME: Currently fixing g to 0
      double f = 1.; // FIXME: Currently fixing f to 1
      vec<double> fGradient;
      fGradient.zeros(path.nD);
      double fLaplacian = 0.;
      //double f = 1. + 2*ZA*(mag_ri_RA);
      //vec<double> fGradient = 2*ZA*((ri_RA/mag_ri_RA));
      //double fLaplacian = 2*ZA*(path.nD-1)*((1./mag_ri_RA));

      // Sum over actions for ri
      std::vector<pair<int,int> > only_ri;
      only_ri.push_back(particles[1]);
      vec<double> actionGradient;
      actionGradient.zeros(path.nD);
      double actionLaplacian = 0.;
      for (auto& action: actionList) {
        actionGradient += action->GetActionGradient(iB,iB+1,only_ri,0);
        actionLaplacian += action->GetActionLaplacian(iB,iB+1,only_ri,0);
      }

      // Sum total
      total += ((g - (1./mag_ri_RA))/(4.*M_PI))*(fLaplacian + f*(-actionLaplacian + dot(actionGradient,actionGradient)) - 2.*dot(fGradient,actionGradient));
    }
  }

  nMeasure += 1;
}

void ContactProbability::Write()
{
  if (nMeasure > 0) {
    // Normalize
    int NA = path.speciesList[iSpeciesA]->nPart;
    int NB = path.speciesList[iSpeciesB]->nPart;
    double norm;
    if (iSpeciesA == iSpeciesB)
      norm = 0.5*nMeasure*NA*(NA-1.)*path.nBead;
    else
      norm = nMeasure*NA*NB*path.nBead;
    total /= norm;

    // Write to file
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet("/"+prefix, "x", total);
    } else {
      out.AppendDataSet("/"+prefix, "x", total);
    }

    Reset();
  }
}
