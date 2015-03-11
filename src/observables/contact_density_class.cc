#include "contact_density_class.h"

void ContactDensity::Init(Input &in)
{
  // Read in species info
  species_A = in.GetAttribute<std::string>("species_A");
  species_B = in.GetAttribute<std::string>("species_B");
  path.GetSpeciesInfo(species_A, species_A_i);
  path.GetSpeciesInfo(species_B, species_B_i);

  // Write things to file
  out.Write(prefix+"/species_A", species_A);
  out.Write(prefix+"/species_B", species_B);

  // Read in Z_A
  Z_A = in.GetAttribute<uint>("Z_A");

  // Generate action list
  std::vector<std::string> species_list;
  species_list.push_back(species_A);
  species_list.push_back(species_B);
  for (auto& action: full_action_list) {
    for (auto& sA: species_list) {
      if (std::find(action->species_list.begin(), action->species_list.end(), sA) != action->species_list.end()) {
        action_list.push_back(action);
        break;
      }
    }
  }

  Reset();
}

void ContactDensity::Reset()
{
  total = 0;
  n_measure = 0;
}

// Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
void ContactDensity::Accumulate()
{
  path.SetMode(1);

  // Form particle pairs
  std::vector<std::vector<std::pair<uint,uint>>> particle_pairs;
  if (species_A_i == species_B_i) { // Homogeneous
    for (uint p_i=0; p_i<path.species_list[species_A_i]->n_part-1; ++p_i) {
      for (uint p_j=p_i+1; p_j<path.species_list[species_B_i]->n_part; ++p_j) {
        std::vector<std::pair<uint,uint>> particles;
        particles.push_back(std::make_pair(species_A_i,p_i));
        particles.push_back(std::make_pair(species_B_i,p_j));
        particle_pairs.push_back(particles);
      }
    }
  } else { // Homologous
    for (uint p_i=0; p_i<path.species_list[species_A_i]->n_part; ++p_i) {
      for (uint p_j=0; p_j<path.species_list[species_B_i]->n_part; ++p_j) {
        std::vector<std::pair<uint,uint>> particles;
        particles.push_back(std::make_pair(species_A_i,p_i));
        particles.push_back(std::make_pair(species_B_i,p_j));
        particle_pairs.push_back(particles);
      }
    }
  }

  // Add up contact probability
  // FIXME: Currently only looking at origin
  for (auto& particles: particle_pairs) {
    for (uint b_i=0; b_i<path.n_bead; ++b_i) {
      // Set r's
      vec<double> RA = path(particles[0].first,particles[0].second,b_i)->r;
      vec<double> ri = path(particles[1].first,particles[1].second,b_i)->r;

      // Get differences
      vec<double> ri_RA(path.Dr(ri, RA));
      double mag_ri_RA = mag(ri_RA);

      // Compute functions
      double g = 0.; // FIXME: Currently fixing g to 0
      double f = 1.; // FIXME: Currently fixing f to 1
      vec<double> gradient_f;
      gradient_f.zeros(path.n_d);
      double laplacian_f = 0.;
      //double f = 1. + 2*Z_A*(mag_ri_RA);
      //vec<double> gradient_f = 2*Z_A*((ri_RA/mag_ri_RA));
      //double laplacian_f = 2*Z_A*(path.n_d-1)*((1./mag_ri_RA));

      // Sum over actions for ri
      std::vector<std::pair<uint,uint>> only_ri;
      only_ri.push_back(particles[1]);
      vec<double> gradient_action;
      gradient_action.zeros(path.n_d);
      double laplacian_action = 0.;
      for (auto& action: action_list) {
        gradient_action += action->GetActionGradient(b_i,b_i+1,only_ri,0);
        laplacian_action += action->GetActionLaplacian(b_i,b_i+1,only_ri,0);
      }

      // Sum total
      total += ((g - (1./mag_ri_RA))/(4.*M_PI))*(laplacian_f + f*(-laplacian_action + dot(gradient_action,gradient_action)) - 2.*dot(gradient_f,gradient_action));
    }
  }

  n_measure += 1;
}

void ContactDensity::Write()
{
  if (n_measure > 0) {
    // Normalize
    uint N_A = path.species_list[species_A_i]->n_part;
    uint N_B = path.species_list[species_B_i]->n_part;
    double norm;
    if (species_A_i == species_B_i)
      norm = 0.5*n_measure*N_A*(N_A-1.)*path.n_bead/path.vol;
    else
      norm = n_measure*N_A*N_B*path.n_bead/path.vol;
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
