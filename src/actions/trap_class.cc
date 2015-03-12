#include "trap_class.h"

void Trap::Init(Input &in)
{
  n_images = in.GetAttribute<int>("n_images");
  omega = in.GetAttribute<double>("omega");
  species = in.GetAttribute<std::string>("species");
  species_list.push_back(species);
  path.GetSpeciesInfo(species,species_i);

  out.Write("/Actions/"+name+"/n_images", n_images);
  out.Write("/Actions/"+name+"/omega", omega);
  out.Write("/Actions/"+name+"/species_i", species_i);
}

double Trap::DActionDBeta()
{
  double tot = 0.;

  for (uint p_i=0; p_i<path.n_part; p_i+=1) {
    for (uint b_i=0; b_i<path.n_bead; b_i+=1) {
      tot += dot(path(species_i,p_i,b_i)->r, path(species_i,p_i,b_i)->r);
    }
  }

  return 0.5*omega*omega*(1. + 3.*path.tau*path.tau*omega*omega/12.)*tot;
}

double Trap::GetAction(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level)
{
  if (level > max_level)
    return 0.;

  bool check = false;
  for (uint p=0; p<particles.size(); ++p) {
    if (particles[p].first == species_i) {
      check = true;
      break;
    }
  }
  if (!check)
    return 0.;

  uint skip = 1<<level;
  double level_tau = skip*path.tau;
  double tot = 0.;
  for (uint p=0; p<particles.size(); ++p) {
    uint s_i = particles[p].first;
    uint p_i = particles[p].second;
    for (uint b_i=b0; b_i<b1; b_i+=skip) {
      vec<double> dr(path.n_d);
      if(path.mode)
        dr = path(s_i,p_i,b_i)->r;
      else
        dr = path(s_i,p_i,b_i)->r_c;
      tot += dot(dr, dr);
    }
  }

  return 0.5*level_tau*omega*omega*(1. + path.tau*path.tau*omega*omega/12.)*tot;
}

void Trap::Write()
{

}
