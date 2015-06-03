#ifndef SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_

#include "../action_class.h"

/// Action class for actions affecting only a single species of particles
class SingleAction : public Action
{
protected:
  double i_4_lambda_tau; ///< 1/(4\lambda\tau)
  double lambda; ///< h_bar^2/2m
  uint32_t n_part; ///< Number of particles in the species affected by the action
  uint32_t species_i; ///< Index of species affected by the action
  std::string species; ///< Name of species affected by the action
public:
  /// Constructor only instatiates parent Action class
  SingleAction(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {
    // Set species variables
    species = in.GetAttribute<std::string>("species");
    out.Write("Actions/"+name+"/species", species);
    species_list.push_back(species);
    path.GetSpeciesInfo(species,species_i);
    n_part = path.species_list[species_i]->n_part;
    lambda = path.species_list[species_i]->lambda;
    i_4_lambda_tau = 1./(4.*lambda*path.tau);
  }
};

#endif // SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_
