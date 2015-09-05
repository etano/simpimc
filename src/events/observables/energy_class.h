#ifndef SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
#define SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_

#include "observable_class.h"

/// Measure the energy of each action
class Energy : public Observable
{
private:
  std::vector<bool> first_sector; ///< Whether or not this is the first sector being written to file
  bool measure_per_sector; ///< Whether or not to measure per permutation section
  bool measure_potential; ///< Whether or not to measure the potential
  bool use_thermal_estimator; ///< Whether or not to use the thermal estimator of the energy
  bool use_virial_estimator; ///< Whether or not to use the virial estimator of the energy
  uint32_t b_i_last; ///< Which bead to average over for virial estimator
  uint32_t virial_window_size; ///< Size of averaging window for virial estimator
  std::vector<std::vector<std::pair<uint32_t,double>>> sector_energies; ///< Vector of vectors of (sector,energy_per_sector) pairs for each species
  std::vector<std::shared_ptr<Action>> &action_list; ///< Reference to vector of all actions
  vec<double> energies; ///< Vector of energies for each action
  vec<double> potentials; ///< Vector of potentials for each action

  /// Thermal estimator of the energy
  double ThermalEstimator()
  {
    double total_energy = 0.;
    int sign = path.CalcSign();
    for (uint32_t i=0; i<action_list.size(); ++i) {
      if (!action_list[i]->is_importance_weight) {
        double action_energy = sign*path.importance_weight*action_list[i]->DActionDBeta();
        energies(i) += action_energy;
        total_energy += action_energy;
      }
    }

    return total_energy;
  }

  /// Virial estimator of the energy
  double VirialEstimator()
  {
    // Kinetic term
    int sign = path.CalcSign();
    double kinetic_energy = 0.;
    for (auto &action: action_list) {
      if (action->type == "Kinetic") {
        kinetic_energy += action->VirialEnergy(virial_window_size);
      }
    }

    // Centroid term
    for (auto &species: path.species_list) {
      if (species->lambda > 0.) {
        uint32_t s_i = species->s_i;

        //// Compute centroids
        //field<vec<double>> centroids(species->n_part,path.n_bead);
        //for (uint32_t p_i=0; p_i<species->n_part; p_i++) {

        //  // First centroid
        //  centroids(p_i,0).zeros(path.n_d);
        //  std::shared_ptr<Bead> b_0 = path(s_i,p_i,0);
        //  std::shared_ptr<Bead> b_1 = path.GetNextBead(path(s_i,p_i,0),1);
        //  for (uint32_t skip=1; skip<virial_window_size; skip++) {
        //    centroids(p_i,0) += path.GetR(b_1);
        //    b_1 = path.GetNextBead(b_1,1);
        //  }

        //  // Other centroids
        //  //
        //  // |_____.___o_|_____.__o__|_o___._____|
        //  //
        //  for (uint32_t b_i=1; b_i<path.n_bead; b_i++) {
        //    centroids(p_i,b_i) = centroids(p_i,b_i-1);
        //    centroids(p_i,b_i) -= path.GetR(b_0);
        //    centroids(p_i,b_i) += path.GetR(b_1);
        //    b_0 = path.GetNextBead(b_0,1);
        //    b_1 = path.GetNextBead(b_1,1);
        //  }

        //  // Normalize
        //  for (uint32_t b_i=0; b_i<path.n_bead; b_i++)
        //    centroids(p_i,b_i) /= virial_window_size; // FIXME: Does this need to be put in the box?
        //}

        // Multiply by force and sum
        b_i_last = path.bead_loop(b_i_last+1);
        double centroid_term = 0.;
        for (uint32_t p_i=0; p_i<species->n_part; p_i++) {
          std::vector<std::pair<uint32_t,uint32_t>> particles;
          particles.push_back(std::make_pair(s_i,p_i));
          for (uint32_t b_i=0; b_i<path.n_bead; b_i++) {
            vec<double> dr(path.Dr(path(s_i,p_i,b_i),path(s_i,p_i,0)));//centroids(p_i,b_i)));
            for (uint32_t i=0; i<action_list.size(); ++i) {
              if (!action_list[i]->is_importance_weight && action_list[i]->type!="Kinetic") {
                vec<double> action_gradient = action_list[i]->GetActionGradient(b_i, b_i+1, particles, 0);
                centroid_term += dot(action_gradient,dr);
              }
            }
          }
        }
        kinetic_energy += centroid_term/(2.*path.tau);

        // Record kinetic energy
        for (uint32_t i=0; i<action_list.size(); ++i) {
          if (action_list[i]->type == "Kinetic") {
            for (auto &species_name : action_list[i]->species_list) {
              if (species_name == species->name)
                 energies(i) += sign*path.importance_weight*kinetic_energy;
            }
          }
        }
      }

    }

    // Interacting term
    double interaction_energy = 0.;
    for (uint32_t i=0; i<action_list.size(); ++i) {
      if (!action_list[i]->is_importance_weight && action_list[i]->type!="Kinetic") {
        double action_energy = sign*path.importance_weight*action_list[i]->DActionDBeta();
        energies(i) += action_energy;
        interaction_energy += action_energy;
      }
    }

    return kinetic_energy + interaction_energy;
  }

  /// Potential estimator
  double PotentialEstimator()
  {
    double total_potential = 0.;
    int sign = path.CalcSign();
    for (uint32_t i=0; i<action_list.size(); ++i) {
      if (!action_list[i]->is_importance_weight) {
        double action_potential = sign*path.importance_weight*action_list[i]->Potential();
        potentials(i) += action_potential;
        total_potential += action_potential;
      }
    }

    return total_potential;
  }

  /// Accumulate the observable
  virtual void Accumulate()
  {
    // Measure energy
    path.SetMode(NEW_MODE);
    double total_energy = 0.;
    if (use_thermal_estimator)
      total_energy = ThermalEstimator();
    else if (use_virial_estimator)
      total_energy = VirialEstimator();

    double total_potential = 0.;
    if (measure_potential)
      total_potential += PotentialEstimator();

    // Store sector with total energy
    if (measure_per_sector) {
      // Loop through species
      for (auto &species: path.species_list) {
        std::pair<uint32_t,double> sector_energy(species->GetPermSector(),species->CalcSign()*total_energy/path.n_bead);
        sector_energies[species->s_i].push_back(sector_energy);
      }
    }

    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_measure = 0;
    energies.zeros();
    if (measure_potential)
      potentials.zeros();

    // Loop through species
    if (measure_per_sector) {
      for (auto &species: path.species_list)
        sector_energies[species->s_i].clear();
    }
  }
public:
  /// Constructor calls Init
  Energy(Path &path, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &out)
    : action_list(tmp_action_list), Observable(path, in, out)
  {
    b_i_last = 0;

    // Whether or not to measure things
    measure_potential = in.GetAttribute<bool>("measure_potential",0);
    use_virial_estimator = in.GetAttribute<bool>("use_virial_estimator",0);
    if (use_virial_estimator) {
      virial_window_size = in.GetAttribute<uint32_t>("virial_window_size",path.n_bead);
      use_thermal_estimator = in.GetAttribute<bool>("use_thermal_estimator",0);
    } else
      use_thermal_estimator = in.GetAttribute<bool>("use_thermal_estimator",1);

    // Initialize permutation sector things if tracking them
    measure_per_sector = in.GetAttribute<uint32_t>("measure_per_sector",0);
    if (measure_per_sector) {
      uint32_t sector_max = in.GetAttribute<uint32_t>("sector_max",0);
      out.CreateGroup(prefix+"sectors");

      // Loop through all species
      sector_energies.resize(path.n_species);
      first_sector.resize(path.n_species);
      for (auto &species: path.species_list) {
        // Set up permutation sectors
        species->SetupPermSectors(sector_max);
        first_sector[species->s_i] = true;

        // Write out possible sectors
        mat<uint32_t> tmp_perms(zeros<mat<uint32_t>>(species->n_part,species->poss_perms.size()));
        std::map<std::vector<uint32_t>,uint32_t>::iterator tmp_iterator;
        for(tmp_iterator = species->poss_perms.begin(); tmp_iterator != species->poss_perms.end(); tmp_iterator++) {
          std::vector<uint32_t> tmpPerm = (*tmp_iterator).first;
          for (uint32_t j=0; j<tmpPerm.size(); ++j)
            tmp_perms(tmpPerm[j]-1,(*tmp_iterator).second)++;
        }
        out.CreateGroup(prefix+"sectors/"+species->name);
        std::string data_type = "avg_pairs";
        out.Write(prefix+"sectors/"+species->name+"/data_type",data_type);
        vec<uint32_t> tmp_perm_indices(species->poss_perms.size());
        for (uint32_t i=0; i<species->poss_perms.size(); ++i)
          tmp_perm_indices(i) = i;
        out.Write(prefix+"sectors/"+species->name+"/x", tmp_perm_indices);
        out.Write(prefix+"sectors/"+species->name+"/poss_perms", tmp_perms);
      }
    }

    // Resize things
    energies.set_size(action_list.size());
    if (measure_potential)
      potentials.set_size(action_list.size());
    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      double norm = path.n_bead*n_measure;

      // Write energies
      energies = energies/norm;
      double E = sum(energies);
      if (first_time) {
        out.CreateGroup(prefix+"total");
        out.CreateExtendableDataSet("/"+prefix+"total/", "x", E);
        std::string data_type = "scalar";
        out.Write(prefix+"total/data_type",data_type);
        for (uint32_t i=0; i<action_list.size(); ++i) {
          out.CreateGroup(prefix+action_list[i]->name);
          out.CreateExtendableDataSet("/"+prefix+action_list[i]->name+"/", "x", energies(i));
          out.Write(prefix+action_list[i]->name+"/data_type", data_type);
        }
      } else {
        out.AppendDataSet("/"+prefix+"total/", "x", E);
        for (uint32_t i=0; i<action_list.size(); ++i)
          out.AppendDataSet("/"+prefix+action_list[i]->name+"/", "x", energies(i));
      }

      // Write potentials
      if (measure_potential) {
        potentials = potentials/norm;
        double V = sum(potentials);
        if (first_time) {
          out.CreateGroup(prefix+"v_total");
          out.CreateExtendableDataSet("/"+prefix+"v_total/", "x", V);
          std::string data_type = "scalar";
          out.Write(prefix+"v_total/data_type",data_type);
          for (uint32_t i=0; i<action_list.size(); ++i) {
            out.CreateGroup(prefix+"v_"+action_list[i]->name);
            out.CreateExtendableDataSet("/"+prefix+"v_"+action_list[i]->name+"/", "x", potentials(i));
            out.Write(prefix+"v_"+action_list[i]->name+"/data_type", data_type);
          }
        } else {
          out.AppendDataSet("/"+prefix+"v_total/", "x", V);
          for (uint32_t i=0; i<action_list.size(); ++i)
            out.AppendDataSet("/"+prefix+"v_"+action_list[i]->name+"/", "x", potentials(i));
        }
      }

      // Write sector energies
      if (measure_per_sector) {
        // Loop through species
        for (auto &species: path.species_list) {
          // Map out the sectors std::vector
          std::map<uint32_t,std::vector<double>> sector_map;
          for (uint32_t i=0; i<n_measure; i++) {
            std::pair<uint32_t,double> sector_energy = sector_energies[species->s_i].back();
            uint32_t sector = sector_energy.first;
            double energy = sector_energy.second;
            if (sector_map.find(sector) == sector_map.end()) {
              std::vector<double> sector_info = {energy,0.,1};
              sector_map.insert(std::pair<uint32_t,std::vector<double> >(sector, sector_info));
            } else {
              double E0 = sector_map[sector][0];
              double var0 = sector_map[sector][1];
              double N0 = sector_map[sector][2];
              double N1 = N0 + 1;
              double E1 = (N0*E0 + energy)/N1; // New average
              double var1 = ((N0*var0 + N0*E0*E0 + energy*energy)/N1) - E1*E1;
              sector_map[sector][0] = E1;
              sector_map[sector][1] = var1;
              sector_map[sector][2] = N1;
            }
            sector_energies[species->s_i].pop_back();
          }

          // Put thestd::map into an array and write
          std::map<uint32_t,uint32_t>::iterator it;
          for(auto& sector_info: sector_map) {
            vec<double> sector_info_vec(4);
            sector_info_vec(0) = sector_info.first;
            sector_info_vec(1) = sector_info.second[0];
            sector_info_vec(2) = sector_info.second[1];
            sector_info_vec(3) = sector_info.second[2];
            if (first_time && first_sector[species->s_i]) {
              first_sector[species->s_i] = false;
              out.CreateExtendableDataSet("/"+prefix+"/sectors/"+species->name+"/", "y", sector_info_vec);
            } else
              out.AppendDataSet("/"+prefix+"/sectors/"+species->name+"/", "y", sector_info_vec);
          }
        }
      }

      if (first_time)
        first_time = 0;

      Reset();
    }
  }

};

#endif // SIMPIMC_OBSERVABLES_ENERGY_CLASS_H_
