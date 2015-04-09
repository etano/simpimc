#ifndef SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_
#define SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_

#include "perm_bisect_class.h"

/// Permuting bisection move that builds the permutation iteratively
class PermBisectIterative : public PermBisect
{
private:
  /// Update the permutation table
  void UpdatePermTable();

  /// Select which cycle to use in the permutation
  bool SelectCycleIterative(Cycle &c);

  /// Attempt the move
  virtual bool Attempt();
public:
  /// Constructor instantiates parent class and calls Init
  PermBisectIterative(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : PermBisect(path, rng, action_list, in, out)
  {}

};

#endif // SIMPIMC_MOVES_PERM_BISECT_ITERATIVE_CLASS_H_
