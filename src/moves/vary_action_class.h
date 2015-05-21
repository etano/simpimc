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
  virtual void Accept();

  /// Attempt the move
  virtual bool Attempt();

  /// Initializes the move
  virtual void Init(Input &in) {};

  /// Rejects the move
  virtual void Reject();

  /// Resets the relevant counters
  virtual void Reset();

public:
  /// Constructor instantiates parent class and calls Init
  VaryAction(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Move(path, rng, action_list, in, out)
  {
    // Get species
    std::string species = in.GetAttribute<std::string>("species");
    path.GetSpeciesInfo(species,species_i);

    // Get relevant varied action
    std::string action_name = in.GetAttribute<std::string>("action_name");

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
  virtual void Write();

};

#endif // SIMPIMC_MOVES_VARY_ACTION_CLASS_H_
