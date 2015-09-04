#ifndef SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
#define SIMPIMC_MOVES_VARY_ACTION_CLASS_H_

#include "move_class.h"

/// Move to vary the parameters of a given action
class VaryAction : public Move
{
private:
  bool switch_param_sets; ///< Whether or not to actually switch parameter sets
  std::shared_ptr<Action> action; ///< Pointer to the variable action
  uint32_t n_param_sets; ///< Number of parameter sets
  uint32_t species_i; ///< Index of a species affected by the action
  uint32_t param_set_new; ///< Index of new parameter set for the action
  uint32_t param_set_old; ///< Index of old parameter set for the action
  mat<double> n_param_accept; ///< Number of parameter switches accepted
  mat<double> n_param_attempt; ///< Number of parameter switches attempted

  /// Accept the move
  virtual void Accept()
  {
    if (switch_param_sets) {
      // Call accept for each action
      action->SetParamSet(param_set_new);
      action->Accept();
    } else {
      // Call accept for each action
      action->SetParamSet(param_set_old);
      action->Reject();
    }

    // Iterate counters
    n_param_accept(param_set_old,param_set_new) += 1;
    n_param_attempt(param_set_old,param_set_new) += 1;
  }

  /// Attempt the move
  virtual bool Attempt()
  {
    // Insert dummy particle
    std::vector<std::pair<uint32_t,uint32_t>> particles;
    particles.push_back(std::make_pair(species_i,0));

    // Get old and set new parameter sets
    param_set_old = action->GetParamSet();
    action->SetParamSet(rng.UnifRand(action->GetNumParamSets())-1);
    param_set_new = action->GetParamSet();

    // Calculate action change
    double old_action = 0.;
    double new_action = 0.;
    if (param_set_old != param_set_new) {
      // Old action
      path.SetMode(OLD_MODE);
      action->SetParamSet(param_set_old);
      old_action += action->GetAction(0, path.n_bead-1, particles, 0);

      // New action
      path.SetMode(NEW_MODE);
      action->SetParamSet(param_set_new);
      new_action += action->GetAction(0, path.n_bead-1, particles, 0);
    }

    double current_action_change = new_action - old_action;
    double log_accept_probablity = -current_action_change;

    // Metropolis reject step
    if (log_accept_probablity < log(rng.UnifRand()))
      return 0;

    return 1;
  }

  /// Rejects the move
  virtual void Reject()
  {
    // Call reject for each action
    action->SetParamSet(param_set_old);
    action->Reject();

    // Iterate counters
    n_param_attempt(param_set_old,param_set_new) += 1;
  }

  /// Resets the relevant counters
  virtual void Reset()
  {
    // Reset counters
    n_param_accept.zeros();
    n_param_attempt.zeros();
  }

public:
  /// Constructor instantiates parent class and reads in relevant inputs
  VaryAction(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    // Get species
    std::string species = in.GetAttribute<std::string>("species");
    path.GetSpeciesInfo(species,species_i);

    // Get relevant varied action
    std::string action_name = in.GetAttribute<std::string>("action_name");

    std::cout << "Setting up " << name << " for " << action_name << " with species " << species << "..." << std::endl;

    // Select action from list
    for (auto& t_action : action_list) {
      if (t_action->name == action_name)
        action = t_action;
    }

    // Whether or not to actually switch parameters
    switch_param_sets = in.GetAttribute<bool>("switch_param_sets",1);

    // Resize things
    n_param_sets = action->GetNumParamSets();
    n_param_accept.set_size(n_param_sets,n_param_sets);
    n_param_attempt.set_size(n_param_sets,n_param_sets);
    Reset();

    // Write entries in output file
    out.CreateGroup("Moves/"+name+"/ParamSets");
    std::string data_type = "matrix";
    out.Write("/Moves/"+name+"/ParamSets/data_type",data_type);
    vec<uint32_t> shape(2);
    shape(0) = n_param_sets;
    shape(1) = n_param_sets;
    out.Write("/Moves/"+name+"/ParamSets/shape",shape);
    out.Write("/Moves/"+name+"/ParamSets/switch_param_sets",switch_param_sets);
  }

  /// Writes relevant information about the move to the output file
  virtual void Write()
  {
    mat<double> param_accept_ratio(n_param_sets,n_param_sets);
    for (uint32_t i=0; i<n_param_sets; i++) {
      for (uint32_t j=0; j<n_param_sets; j++) {
        double ratio = n_param_accept(i,j)/n_param_attempt(i,j);
        param_accept_ratio(i,j) = ratio;
      }
    }

    // Write
    if (first_time) {
      out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "n_param_accept", n_param_accept);
      out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "n_param_attempt", n_param_attempt);
      out.CreateExtendableDataSet("/Moves/"+name+"/ParamSets/", "y", param_accept_ratio);
    } else {
      out.AppendDataSet("/Moves/"+name+"/ParamSets/", "n_param_attempt", n_param_attempt);
      out.AppendDataSet("/Moves/"+name+"/ParamSets/", "n_param_accept", n_param_accept);
      out.AppendDataSet("/Moves/"+name+"/ParamSets/", "y", param_accept_ratio);
    }

    // Reset
    Reset();

    Move::Write();
  }

};

#endif // SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
