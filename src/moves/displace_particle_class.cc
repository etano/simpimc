#include "displace_particle_class.h"

void DisplaceParticle::Init(Input &in)
{
  species = in.GetAttribute<std::string>("species");
  path.GetSpeciesInfo(species,species_i);
  step_size = in.GetAttribute<double>("step_size",path.L/10.);

  // Generate action list
  std::vector<std::string> species_list;
  species_list.push_back(species);
  GenerateActionList(species_list);
}

// Accept current move
void DisplaceParticle::Accept()
{
  // Move Accepted, so copy new coordinates
  path.StoreR(affected_beads);
  path.StoreRhoKP(affected_beads);
  for (uint32_t b_i=0; b_i<path.n_bead; ++b_i)
    path.StoreRhoK(b_i,species_i);

  // Call accept for each action
  for (auto& action: action_list)
    action->Accept();
}

// Reject current move
void DisplaceParticle::Reject()
{
  // Move rejected, so return old coordinates
  path.RestoreR(affected_beads);
  path.RestoreRhoKP(affected_beads);
  for (uint32_t b_i=0; b_i<path.n_bead; ++b_i)
    path.RestoreRhoK(b_i,species_i);

  // Call reject for each action
  for (auto& action: action_list)
    action->Reject();
}

// Displace particle Move
bool DisplaceParticle::Attempt()
{
  // Set which particles are affected by the move
  uint32_t p_i = rng.UnifRand(path.species_list[species_i]->n_part) - 1;  // Pick particle at random
  std::vector<std::pair<uint32_t,uint32_t>> particles;
  std::pair<uint32_t,uint32_t> particle(species_i,p_i);
  particles.push_back(particle);

  // New sampling
  path.SetMode(NEW_MODE);
  vec<double> dr(path.n_d);
  rng.UnifRand(dr, step_size);

  // Set which beads are affected by the move
  // and move them
  affected_beads.clear();
  for(uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
    affected_beads.push_back(path(species_i,p_i,b_i));
    path.GetR(path(species_i,p_i,b_i)) += dr;
  }

  // Calculate action change
  double old_action = 0.;
  double new_action = 0.;
  for (auto& action: action_list) {
    // Old action
    path.SetMode(OLD_MODE);
    old_action += action->GetAction(0, path.n_bead, particles, 0);

    // New action
    path.SetMode(NEW_MODE);
    new_action += action->GetAction(0, path.n_bead, particles, 0);
  }

  double log_accept_probablity = old_action - new_action;

  // Metropolis reject step
  if (log_accept_probablity < log(rng.UnifRand()))
    return 0;
  else
    return 1;
}
