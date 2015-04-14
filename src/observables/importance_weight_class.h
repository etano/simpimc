#ifndef SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_
#define SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_

#include "observable_class.h"

/// Accumulates any and all importance weights
class ImportanceWeight : public Observable
{
private:
  std::vector<std::shared_ptr<Action>> &action_list; /// List of all actions
  vec<double> importance_weights; ///< Vector of importance weights

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  ImportanceWeight(Path &path, std::vector<std::shared_ptr<Action>> &t_action_list, Input &in, IO &out)
    : action_list(t_action_list), Observable(path, in, out)
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_
