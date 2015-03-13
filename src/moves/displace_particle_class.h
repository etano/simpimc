#ifndef SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
#define SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_

#include "move_class.h"

class DisplaceParticle : public Move
{
private:
  std::string species;
  uint32_t species_i;
  std::vector<std::shared_ptr<Bead>> affected_beads;
  double step_size;
protected:

public:
  // Constructor
  DisplaceParticle(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif // SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
