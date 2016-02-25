#ifndef SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_
#define SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_

#include "observable_class.h"

/// Accumulates any and all importance weights
class ImportanceWeight : public Observable {
   private:
    std::vector<std::shared_ptr<Action>> &action_list;  /// List of all actions
    vec<double> importance_weights;                     ///< Vector of importance weights

    /// Accumulate the observable
    virtual void Accumulate() {
        path.SetMode(NEW_MODE);
        double importance_weight = 1.;
        for (uint32_t i = 0; i < action_list.size(); ++i) {
            double action_importance_weight = 1.;
            if (action_list[i]->is_importance_weight)
                action_importance_weight = action_list[i]->ImportanceWeight();
            importance_weights(i) += action_importance_weight;
            importance_weight *= action_importance_weight;
        }
        path.SetImportanceWeight(importance_weight);
        n_measure += 1;
    }

    /// Reset the observable's counters
    virtual void Reset() {
        n_measure = 0;
        importance_weights.zeros();
    }

   public:
    /// Constructor calls Init
    ImportanceWeight(Path &path, std::vector<std::shared_ptr<Action>> &t_action_list, Input &in, IO &out)
        : action_list(t_action_list), Observable(path, in, out) {
        importance_weights.set_size(action_list.size());
        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write() {
        if (n_measure > 0) {
            double norm = n_measure;

            // Write importance_weights
            importance_weights = importance_weights / norm;
            double IW = prod(importance_weights);
            if (first_time) {
                out.CreateGroup(prefix + "Total");
                out.CreateExtendableDataSet("/" + prefix + "Total/", "x", IW);
                std::string data_type = "scalar";
                out.Write(prefix + "Total/data_type", data_type);
                for (uint32_t i = 0; i < action_list.size(); ++i) {
                    out.CreateGroup(prefix + action_list[i]->name);
                    out.CreateExtendableDataSet("/" + prefix + action_list[i]->name + "/", "x", importance_weights(i));
                    out.Write(prefix + action_list[i]->name + "/data_type", data_type);
                }
                first_time = 0;
            } else {
                out.AppendDataSet("/" + prefix + "Total/", "x", IW);
                for (uint32_t i = 0; i < action_list.size(); ++i)
                    out.AppendDataSet("/" + prefix + action_list[i]->name + "/", "x", importance_weights(i));
            }

            Reset();
        }
    }
};

#endif  // SIMPIMC_OBSERVABLES_IMPORTANCE_WEIGHT_CLASS_H_
