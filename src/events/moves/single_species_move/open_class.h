#ifndef SIMPIMC_MOVES_OPEN_CLASS_H_
#define SIMPIMC_MOVES_OPEN_CLASS_H_

#include "single_species_move_class.h"

/// This moves displaces whole particles by a defined step size
class Open : public SingleSpeciesMove
{
private:
  /// Accept the move
  virtual void Accept() {}

  /// Attempt the move
  virtual bool Attempt()
  {
    // Pick particle and bead at random
    uint32_t p_i = rng.UnifRand(path.species_list[species_i]->n_part) - 1;
    uint32_t b_i = rng.UnifRand(path.n_bead) - 1;
    std::cout << p_i << " " << b_i << std::endl;

    // Open path at (p_i,b_i)
    path(species_i,p_i,b_i)->next->prev = nullptr;
    path(species_i,p_i,b_i)->next_c->prev_c = nullptr;
    path(species_i,p_i,b_i)->next = nullptr;
    path(species_i,p_i,b_i)->next_c = nullptr;

    return 1;
  }

  /// Rejects the move
  virtual void Reject() {}
public:
  /// Constructor instantiates parent
  Open(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {}

};

#endif // SIMPIMC_MOVES_OPEN_CLASS_H_
