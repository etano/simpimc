#include "energy_class.h"

void Energy::Init(Input &in)
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

    // Loop through all species
    sector_energies.resize(path.n_species);
    for (auto &species: path.species_list) {
      // Set up permutation sectors
      uint32_t species_i = species->s_i;
      path.SetupPermSectors(species->s_i, sector_max);
      first_sector = true;

      // Write out possible sectors
      mat<uint32_t> tmp_perms;
      tmp_perms.zeros(species->n_part,path.poss_perms[species_i].size());
      std::map<std::vector<uint32_t>,uint32_t>::iterator tmp_iterator;
      for(tmp_iterator = path.poss_perms[species_i].begin(); tmp_iterator != path.poss_perms[species_i].end(); tmp_iterator++) {
        std::vector<uint32_t> tmpPerm = (*tmp_iterator).first;
        for (uint32_t j=0; j<tmpPerm.size(); ++j)
          tmp_perms(tmpPerm[j]-1,(*tmp_iterator).second)++;
      }
      out.CreateGroup(prefix+"sectors");
      out.CreateGroup(prefix+"sectors/"+species->name);
      std::string data_type = "avg_pairs";
      out.Write(prefix+"sectors/"+species->name+"/data_type",data_type);
      vec<uint32_t> tmp_perm_indices(path.poss_perms[species_i].size());
      for (uint32_t i=0; i<path.poss_perms[species_i].size(); ++i)
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

void Energy::Reset()
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

void Energy::Accumulate()
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
      uint32_t species_i = species->s_i;
      std::vector<uint32_t> cycle;
      path.SetCycleCount(species_i, cycle);
      uint32_t sector = path.GetPermSector(species_i, cycle);
      std::pair<uint32_t,double> sector_energy(sector,path.sign*total_energy/path.n_bead);
      sector_energies[species_i].push_back(sector_energy);
    }
  }

  n_measure += 1;
}

double Energy::ThermalEstimator()
{
  double total_energy = 0.;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    if (!action_list[i]->is_importance_weight) {
      double action_energy = path.sign*path.importance_weight*action_list[i]->DActionDBeta();
      energies(i) += action_energy;
      total_energy += action_energy;
    }
  }

  return total_energy;
}

double Energy::PotentialEstimator()
{
  double total_potential = 0.;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    if (!action_list[i]->is_importance_weight) {
      double action_potential = path.sign*path.importance_weight*action_list[i]->Potential();
      potentials(i) += action_potential;
      total_potential += action_potential;
    }
  }

  return total_potential;
}

double Energy::VirialEstimator()
{
  // Kinetic term
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
               energies(i) += path.sign*path.importance_weight*kinetic_energy;
          }
        }
      }
    }

  }

  // Interacting term
  double interaction_energy = 0.;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    if (!action_list[i]->is_importance_weight && action_list[i]->type!="Kinetic") {
      double action_energy = path.sign*path.importance_weight*action_list[i]->DActionDBeta();
      energies(i) += action_energy;
      interaction_energy += action_energy;
    }
  }

  return kinetic_energy + interaction_energy;
}

void Energy::Write()
{
  if (n_measure > 0) {
    double norm = path.n_bead*n_measure;

    // Write energies
    energies = energies/norm;
    double E = sum(energies);
    if (first_time) {
      out.CreateGroup(prefix+"Total");
      out.CreateExtendableDataSet("/"+prefix+"Total/", "x", E);
      std::string data_type = "scalar";
      out.Write(prefix+"Total/data_type",data_type);
      for (uint32_t i=0; i<action_list.size(); ++i) {
        out.CreateGroup(prefix+action_list[i]->name);
        out.CreateExtendableDataSet("/"+prefix+action_list[i]->name+"/", "x", energies(i));
        out.Write(prefix+action_list[i]->name+"/data_type", data_type);
      }
    } else {
      out.AppendDataSet("/"+prefix+"Total/", "x", E);
      for (uint32_t i=0; i<action_list.size(); ++i)
        out.AppendDataSet("/"+prefix+action_list[i]->name+"/", "x", energies(i));
    }

    // Write potentials
    if (measure_potential) {
      potentials = potentials/norm;
      double V = sum(potentials);
      if (first_time) {
        out.CreateGroup(prefix+"VTotal");
        out.CreateExtendableDataSet("/"+prefix+"VTotal/", "x", V);
        std::string data_type = "scalar";
        out.Write(prefix+"VTotal/data_type",data_type);
        for (uint32_t i=0; i<action_list.size(); ++i) {
          out.CreateGroup(prefix+"V"+action_list[i]->name);
          out.CreateExtendableDataSet("/"+prefix+"V"+action_list[i]->name+"/", "x", potentials(i));
          out.Write(prefix+"V"+action_list[i]->name+"/data_type", data_type);
        }
      } else {
        out.AppendDataSet("/"+prefix+"VTotal/", "x", V);
        for (uint32_t i=0; i<action_list.size(); ++i)
          out.AppendDataSet("/"+prefix+"V"+action_list[i]->name+"/", "x", potentials(i));
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
          if (first_time && first_sector) {
            first_sector = false;
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
