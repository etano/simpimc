#ifndef SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_
#define SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_

#include "observable_class.h"

/// Record the fraction of Monte Carlo steps spent with each parameter set of an optimizable action
class RecordOptimizedAction : public Observable
{
private:
  std::shared_ptr<Action> action; ///< Pointer to action being optimized
  std::vector<std::shared_ptr<Action>> &action_list; ///< Reference to list of all actions
  vec<double> param_set_count; ///< Vector of counters for each possible parameter set

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  RecordOptimizedAction(Path &path, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &out)
    : action_list(tmp_action_list), Observable(path, in, out, "histogram")
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_
