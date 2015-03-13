#ifndef SIMPIMC_ACTIONS_TRAP_CLASS_H_
#define SIMPIMC_ACTIONS_TRAP_CLASS_H_

#include "action_class.h"

class Trap : public Action
{
private:
  int n_images;
  uint32_t species_i, max_level;
  double omega;
  std::string species;
protected:

public:
  // Constructor
  Trap(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual double DActionDBeta();
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);
  virtual void Write();
};

#endif // SIMPIMC_ACTIONS_TRAP_CLASS_H_
