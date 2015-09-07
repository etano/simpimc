#ifndef SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_
#define SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_

#include "../move_class.h"

/// Parent class for all moves involving only a single species
class SingleSpeciesMove : public Move
{
protected:
  double i_4_lambda_tau; ///< 1/(4\lambda\tau) for the affected species
  std::shared_ptr<Species> species; ///< Pointer to relevant species
public:
  /// Constructor gets information about the species from the input and generates the action list
  SingleSpeciesMove(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &t_action_list, Input &in, IO &out)
    : Move(path,rng,t_action_list,in,out)
  {
    // Get species
    std::string species_name = in.GetAttribute<std::string>("species");
    species = path.GetSpecies(species_name);
    i_4_lambda_tau = 1./(4.*species->GetLambda()*path.GetTau());

    // Generate action list
    GenerateActionList(t_action_list,species);

    std::cout << "Setting up " << name << " for " << species_name << "..." << std::endl;
  }
};

#endif // SIMPIMC_MOVES_SINGLE_SPECIES_MOVE_CLASS_H_
