#ifndef SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
#define SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_

#include "single_species_move_class.h"

/// This move attempts to shift the reference slice
class ShiftRefSlice : public SingleSpeciesMove
{
private:
  uint32_t ref_bead_0; ///< Original refence slice
  uint32_t ref_bead_1; ///< New reference slice

  /// Accept the move
  virtual void Accept();

  /// Attempt the move
  virtual bool Attempt();

  /// Rejects the move
  virtual void Reject();
public:
  /// Constructor only instantiates parent
  ShiftRefSlice(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {}
};

#endif // SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
