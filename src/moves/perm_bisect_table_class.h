#ifndef SIMPIMC_MOVES_PERM_TABLE_BISECT_CLASS_H_
#define SIMPIMC_MOVES_PERM_TABLE_BISECT_CLASS_H_

#include "perm_bisect_class.h"

/// Permuting bisection move which samples a permutation table
class PermBisectTable : public PermBisect
{
private:
  /// Construct the permutation table
  double ConstructPermTable();

  /// Update the permutation table
  void UpdatePermTable();

  /// Build all possible cycles
  void BuildCycles();

  /// Select which cycle to perform a permutation with
  uint32_t SelectCycle(const double permTot);

  /// Attempt the move
  virtual bool Attempt();
public:
  /// Constructor instantiates parent class, sets the number of particles involved in a permutation, and builds all possible cycles
  PermBisectTable(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : PermBisect(path, rng, action_list, in, out)
  {
    // Read in things
    n_perm_part = in.GetAttribute<uint32_t>("n_perm_part");

    // Generate all cycles
    BuildCycles();
  }
};


#endif // SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
