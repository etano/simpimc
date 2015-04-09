#ifndef SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
#define SIMPIMC_MOVES_VARY_ACTION_CLASS_H_

#include "move_class.h"

/// Move to vary the parameters of a given action
class VaryAction : public Move
{
private:
  std::shared_ptr<Action> action; ///< Pointer to the variable action
  uint32_t species_i; ///< Index of a species affected by the action
  uint32_t param_set_new; ///< Index of new parameter set for the action
  uint32_t param_set_old; ///< Index of old parameter set for the action

  /// Accept the move
  virtual void Accept();

  /// Attempt the move
  virtual bool Attempt();

  /// Initializes the move
  virtual void Init(Input &in);

  /// Rejects the move
  virtual void Reject();
public:
  /// Constructor instantiates parent class and calls Init
  VaryAction(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }
};

#endif // SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
