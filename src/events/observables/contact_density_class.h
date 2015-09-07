#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_

#include "observable_class.h"

/// Measures the contact density between two species of particles. Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
class ContactDensity : public Observable
{
private:
  double total; ///< Running total
  uint32_t z_a; ///< Charge of ion-like particle
  std::shared_ptr<Species> species_a; ///< ion species
  std::shared_ptr<Species> species_b; ///< other species
  std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that involves species_a or species_b
  std::vector<std::shared_ptr<Action>> &full_action_list; ///< Vector of pointers to all actions

  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);

    // Form particle pairs
    std::vector<std::pair<uint32_t,uint32_t>> particle_pairs;
    if (species_a == species_b) { // Homogeneous
      for (uint32_t p_i=0; p_i<species_a->GetNPart()-1; ++p_i)
        for (uint32_t p_j=p_i+1; p_j<species_b->GetNPart(); ++p_j)
          particle_pairs.push_back(std::make_pair(p_i,p_j));
    } else { // Homologous
      for (uint32_t p_i=0; p_i<species_a->GetNPart(); ++p_i)
        for (uint32_t p_j=0; p_j<species_b->GetNPart(); ++p_j)
          particle_pairs.push_back(std::make_pair(p_i,p_j));
    }

    // Add up contact probability
    // TODO: Currently only looking at origin
    double tot = 0.;
    size_t n_particle_pairs(particle_pairs.size());
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t pp_i=0; pp_i<n_particle_pairs; ++pp_i) {
      for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {
        // Set r's
        vec<double> RA = species_a->GetBead(particle_pairs[pp_i].first,b_i)->GetR();
        vec<double> ri = species_b->GetBead(particle_pairs[pp_i].second,b_i)->GetR();

        // Get differences
        vec<double> ri_RA(path.Dr(ri, RA));
        double mag_ri_RA = mag(ri_RA);

        // Compute functions
        double g = 0.; // TODO: Currently fixing g to 0
        double f = 1.; // TODO: Currently fixing f to 1
        vec<double> gradient_f(zeros<vec<double>>(path.GetND()));
        double laplacian_f = 0.;
        //double f = 1. + 2*z_a*(mag_ri_RA);
        //vec<double> gradient_f = 2*z_a*((ri_RA/mag_ri_RA));
        //double laplacian_f = 2*z_a*(path.GetND()-1)*((1./mag_ri_RA));

        // Sum over actions for ri
        std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> only_ri{std::make_pair(species_a,particle_pairs[pp_i].second)};
        vec<double> gradient_action(zeros<vec<double>>(path.GetND()));
        double laplacian_action = 0.;
        for (auto& action: action_list) {
          gradient_action += action->GetActionGradient(b_i,b_i+1,only_ri,0);
          laplacian_action += action->GetActionLaplacian(b_i,b_i+1,only_ri,0);
        }

        // Sum total
        tot += ((g - (1./mag_ri_RA))/(4.*M_PI))*(laplacian_f + f*(-laplacian_action + dot(gradient_action,gradient_action)) - 2.*dot(gradient_f,gradient_action));
      }
    }

    double cofactor = path.GetSign()*path.GetImportanceWeight();
    total += cofactor*tot;
    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    total = 0;
    n_measure = 0;
  }

public:
  /// Constructor calls Init
  ContactDensity(Path &path, std::vector<std::shared_ptr<Action>>& t_action_list, Input &in, IO &out)
    : full_action_list(t_action_list), Observable(path, in, out, "scalar")
  {
    // Read in species info
    std::string species_a_name = in.GetAttribute<std::string>("species_a");
    std::string species_b_name = in.GetAttribute<std::string>("species_b");
    species_a = path.GetSpecies(species_a_name);
    species_b = path.GetSpecies(species_b_name);

    // Write things to file
    out.Write(prefix+"/species_a", species_a_name);
    out.Write(prefix+"/species_b", species_b_name);

    // Read in z_a
    z_a = in.GetAttribute<uint32_t>("z_a");

    // Generate action list
    for (auto& action: full_action_list)
        if ((std::find(action->species_list.begin(), action->species_list.end(), species_a)!=action->species_list.end()) or (std::find(action->species_list.begin(), action->species_list.end(), species_b)!=action->species_list.end()))
          action_list.push_back(action);

    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      // Normalize
      double norm;
      if (species_a == species_b)
        norm = 0.5*n_measure*species_a->GetNPart()*(species_a->GetNPart()-1.)*path.GetNBead()/path.GetVol();
      else
        norm = n_measure*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead()/path.GetVol();
      total /= norm;

      // Write to file
      if (first_time) {
        first_time = 0;
        out.CreateExtendableDataSet("/"+prefix, "x", total);
      } else {
        out.AppendDataSet("/"+prefix, "x", total);
      }

      Reset();
    }
  }
};

#endif // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
