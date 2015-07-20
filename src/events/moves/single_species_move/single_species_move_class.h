#ifndef SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_
#define SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_

#include "../move_class.h"

/// Parent class for all moves involving only a single species
class SingleSpeciesMove : public Move
{
protected:
  double lambda; ///< \hbar^2/2m for the affected species
  double i_4_lambda_tau; ///< 1/(4\lambda\tau) for the affected species
  uint32_t n_part; ///< Number of particles in affected species
  uint32_t species_i; ///< Index of affected species
  std::string species; ///< Name of affected species
public:
  /// Constructor gets information about the species from the input and generates the action list
  SingleSpeciesMove(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &t_action_list, Input &in, IO &out)
    : Move(path,rng,t_action_list,in,out)
  {
    // Get species
    species = in.GetAttribute<std::string>("species");
    path.GetSpeciesInfo(species,species_i);
    n_part = path.species_list[species_i]->n_part;
    lambda = path.species_list[species_i]->lambda;
    i_4_lambda_tau = 1./(4.*lambda*path.tau);

    // Generate action list
    std::vector<std::string> species_list;
    species_list.push_back(species);
    GenerateActionList(t_action_list,species_list);

    std::cout << "Setting up " << name << " for " << species << "..." << std::endl;
  }
};

#endif // SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_
