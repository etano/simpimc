#ifndef SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
#define SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_

#include "single_species_move_class.h"

/// This moves displaces whole particles by a defined step size
class DisplaceParticle : public SingleSpeciesMove
{
private:
  double step_size; ///< Distanced displaced during move
  std::vector<std::shared_ptr<Bead>> affected_beads; ///< Beads affected by the move

  /// Accept the move
  virtual void Accept()
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

  /// Attempt the move
  virtual bool Attempt()
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
    std::shared_ptr<Bead> beadA(path(species_i,p_i,0));
    std::shared_ptr<Bead> beadF(path.GetNextBead(beadA,path.n_bead));
    while(beadA != beadF) {
      affected_beads.push_back(beadA);
      path.GetR(beadA) += dr;
      std::shared_ptr<Bead> beadA(path.GetNextBead(beadA,1));
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

  /// Rejects the move
  virtual void Reject()
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
public:
  /// Constructor instantiates parent and sets the step size
  DisplaceParticle(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {
    step_size = in.GetAttribute<double>("step_size",path.L/10.);
  }

};

#endif // SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
