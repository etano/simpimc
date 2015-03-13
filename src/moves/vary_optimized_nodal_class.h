#ifndef SIMPIMC_MOVES_VARY_OPTIMIZED_NODAL_CLASS_H_
#define SIMPIMC_MOVES_VARY_OPTIMIZED_NODAL_CLASS_H_

#include "move_class.h"

class VaryOptimizedNodal : public Move
{
private:
  uint32_t species_i, param_set_0, param_set_1;
protected:

public:
  // Constructor
  VaryOptimizedNodal(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif // SIMPIMC_MOVES_VARY_OPTIMIZED_NODAL_CLASS_H_
