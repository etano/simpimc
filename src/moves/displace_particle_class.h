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
  virtual void Accept();

  /// Attempt the move
  virtual bool Attempt();

  /// Rejects the move
  virtual void Reject();
public:
  /// Constructor instantiates parent and sets the step size
  DisplaceParticle(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {
    step_size = in.GetAttribute<double>("step_size",path.L/10.);
  }

};

#endif // SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
