#ifndef SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
#define SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_

#include "move_class.h"

class ShiftRefSlice : public Move
{
private:
  std::string species;
  uint32_t species_i, ref_bead_0, ref_bead_1;
protected:

public:
  // Constructor
  ShiftRefSlice(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();

};

#endif // SIMPIMC_MOVES_SHIFT_REF_SLICE_CLASS_H_
