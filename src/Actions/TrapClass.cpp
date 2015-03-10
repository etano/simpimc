#include "TrapClass.h"

void Trap::Init(Input &in)
{
  nImages = in.getAttribute<int>("nImages");
  omega = in.getAttribute<double>("omega");
  species = in.getAttribute<string>("species");
  speciesList.push_back(species);
  path.GetSpeciesInfo(species,iSpecies);

  out.Write("/Actions/"+name+"/nImages", nImages);
  out.Write("/Actions/"+name+"/omega", omega);
  out.Write("/Actions/"+name+"/iSpecies", iSpecies);
}

double Trap::DActionDBeta()
{
  double tot = 0.;

  for (uint iP=0; iP<path.nPart; iP+=1) {
    for (uint iB=0; iB<path.nBead; iB+=1) {
      tot += dot(path(iSpecies,iP,iB)->r, path(iSpecies,iP,iB)->r);
    }
  }

  return 0.5*omega*omega*(1. + 3.*path.tau*path.tau*omega*omega/12.)*tot;
}

double Trap::GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> > &particles, const uint level)
{
  if (level > maxLevel)
    return 0.;

  bool check = false;
  for (uint p=0; p<particles.size(); ++p) {
    if (particles[p].first == iSpecies) {
      check = true;
      break;
    }
  }
  if (!check)
    return 0.;

  uint skip = 1<<level;
  double levelTau = skip*path.tau;
  double tot = 0.;
  for (uint p=0; p<particles.size(); ++p) {
    uint iS = particles[p].first;
    uint iP = particles[p].second;
    for (uint iB=b0; iB<b1; iB+=skip) {
      vec<double> dr(path.nD);
      if(path.mode)
        dr = path(iS,iP,iB)->r;
      else
        dr = path(iS,iP,iB)->rC;
      tot += dot(dr, dr);
    }
  }

  return 0.5*levelTau*omega*omega*(1. + path.tau*path.tau*omega*omega/12.)*tot;
}

void Trap::Write()
{

}
