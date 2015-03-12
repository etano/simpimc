#ifndef SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
#define SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_

#include "observable_class.h"

class Energy : public Observable
{
private:
  std::vector<std::shared_ptr<Action>> &action_list;
  vec<double> energies, potentials;
  std::vector<std::pair<uint,double>> sector_energies;
  bool measure_potential, measure_per_sector, first_sector;
  uint species_i;
protected:
public:
  Energy(Path &path, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &out)
    : action_list(tmp_action_list), Observable(path, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
