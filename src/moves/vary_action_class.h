#ifndef SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
#define SIMPIMC_MOVES_VARY_ACTION_CLASS_H_

#include "move_class.h"

class VaryAction : public Move
{
private:
  uint32_t species_i, param_set_old, param_set_new;
protected:

public:
  // Constructor
  VaryAction(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif // SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
