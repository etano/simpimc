#ifndef SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_
#define SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_

#include "observable_class.h"

/// Record the fraction of Monte Carlo steps spent with each parameter set of an optimizable action
class RecordOptimizedAction : public Observable {
   private:
    std::shared_ptr<Action> action;                     ///< Pointer to action being optimized
    std::vector<std::shared_ptr<Action>> &action_list;  ///< Reference to list of all actions
    vec<double> param_set_count;                        ///< Vector of counters for each possible parameter set

    /// Accumulate the observable
    virtual void Accumulate() {
        path.SetMode(NEW_MODE);
        param_set_count(action->GetParamSet()) += 1.;
        n_measure++;
    }

    /// Reset the observable's counters
    virtual void Reset() {
        n_measure = 0;
        param_set_count.zeros();
    }

   public:
    /// Constructor calls Init
    RecordOptimizedAction(Path &path, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &out)
        : action_list(tmp_action_list), Observable(path, in, out, "histogram") {
        // Read in action name
        std::string action_name = in.GetAttribute<std::string>("action_name");

        // Select action from list
        for (auto &t_action : action_list) {
            if (t_action->name == action_name)
                action = t_action;
        }

        // Set up param_set_count
        param_set_count.set_size(action->GetNumParamSets());
        vec<uint32_t> param_set_indices(param_set_count.size());
        for (uint32_t i = 0; i < param_set_count.size(); ++i)
            param_set_indices(i) = i;
        out.Write(prefix + "x", param_set_indices);

        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write() {
        if (n_measure > 0) {
            double norm = 1. * n_measure;

            // Write sign
            param_set_count /= norm;
            if (first_time) {
                first_time = 0;
                out.CreateExtendableDataSet(prefix, "y", param_set_count);
            } else {
                out.AppendDataSet(prefix, "y", param_set_count);
            }

            Reset();
        }
    }
};

#endif  // SIMPIMC_OBSERVABLES_RECORD_OPTIMIZED_ACTION_CLASS_H_
