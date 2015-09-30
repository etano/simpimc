#ifndef SIMPIMC_ACTIONS_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_ACTION_CLASS_H_

#include "../data_structures/path_class.h"

/// Parent class for all actions
class Action
{
protected:
  int n_images; ///< Number of images to include in calculation
  uint32_t param_set_i; ///< Current parameter set being used for model_i
  uint32_t max_level; ///< Maximum bisection level affected by the action
  uint32_t model_i; ///< Current model being used
  std::vector<std::vector<double>> param_sets; ///< Vector of variational
  IO &out; ///< Reference to ouput file
  Path& path; ///< Reference to path
public:
  bool is_importance_weight; ///< Whether or not the action is actually an importance weight
  std::string name; ///< Action name
  std::string type; ///< Action type
  std::vector<std::shared_ptr<Species>> species_list; ///< List of species affected by the action

  /// Constructor sets the name and type, as well as the relevant references to the path and output file. It also creates the space in the output file for which to write information about the action
  Action(Path &tmp_path, Input &in, IO &tmpOut)
    : path(tmp_path), out(tmpOut), is_importance_weight(0)
  {
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    out.CreateGroup("Actions/"+name);
    out.Write("Actions/"+name+"/type",type);

    n_images = in.GetAttribute<int>("n_images",0);
    max_level = in.GetAttribute<uint32_t>("max_level",0);
    out.Write("Actions/"+name+"/n_images",n_images);
    out.Write("Actions/"+name+"/max_level", max_level);
  }

  /// Returns the beta derivative of the action for the whole path
  virtual double DActionDBeta() { return 0.; };

  /// Returns the virial contribution of the action
  virtual double VirialEnergy(const double virial_window_size) { return DActionDBeta(); }

  /// Returns the value of the action between time slices b0 and b1 for a vector of particles
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level) { return 0.; };

  /// Returns the spatial gradient of the action between time slices b0 and b1 for a vector of particles
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level) { return zeros<vec<double>>(path.GetND()); };

  /// Returns the spatial laplacian of the action between time slices b0 and b1 for a vector of particles
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> &particles, const uint32_t level) { return 0.; };

  /// Returns the potential of the action for the whole path
  virtual double Potential() { return 0.; };

  /// Returns the importance weight of the action for the whole path
  virtual double ImportanceWeight() { return 0; };

  /// Write information about the action
  virtual void Write() {};

  /// Accepts relevant information about the action
  virtual void Accept() {};

  /// Rejects relevant information about the action
  virtual void Reject() {};

  /// Returns current parameter set for action
  virtual uint32_t GetParamSet() { return param_set_i; };

  /// Returns the number of parameter sets for the action
  virtual uint32_t GetNumParamSets() { return param_sets.size(); };

  /// Sets the parameter set of the action
  virtual void SetParamSet(uint32_t t_param_set_i) { param_set_i = t_param_set_i; };

};

#endif // SIMPIMC_ACTIONS_ACTION_CLASS_H_
