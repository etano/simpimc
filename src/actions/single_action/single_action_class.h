#ifndef SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_

#include "../action_class.h"

/// Action class for actions affecting only a single species of particles
class SingleAction : public Action
{
protected:
  double i_4_lambda_tau; ///< 1/(4\lambda\tau)
  std::shared_ptr<Species> species; ///< Pointer to relevant species
public:
  /// Constructor only instatiates parent Action class
  SingleAction(Path &path, Input &in, IO &out)
    : Action(path,in,out)
  {
    // Get species
    std::string species_name = in.GetAttribute<std::string>("species");
    species = path.GetSpecies(species_name);
    i_4_lambda_tau = 1./(4.*species->GetLambda()*path.GetTau());
    species_list.push_back(species);
    out.Write("Actions/"+name+"/species", species_name);

    std::cout << "Setting up " << name << " for " << species_name << "..." << std::endl;
  }
};

#endif // SIMPIMC_ACTIONS_SINGLE_ACTION_CLASS_H_
