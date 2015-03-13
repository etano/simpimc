#include "energy_class.h"

void Energy::Init(Input &in)
{
  measure_potential = in.GetAttribute<uint32_t>("measure_potential",0);

  // Initialize permutation sector things if tracking them
  measure_per_sector = in.GetAttribute<uint32_t>("measure_per_sector",0);
  if (measure_per_sector) {
    uint32_t sector_max = in.GetAttribute<uint32_t>("sector_max",0);
    std::string species = in.GetAttribute<std::string>("species");

    // Set up permutation sectors
    path.GetSpeciesInfo(species, species_i);
    path.SetupPermSectors(path.species_list[species_i]->n_part, sector_max);
    first_sector = true;

    // Write out possible sectors
    mat<uint32_t> tmp_perms;
    tmp_perms.zeros(path.species_list[species_i]->n_part,path.poss_perms.size());
    std::map<std::vector<uint32_t>,uint32_t>::iterator tmp_iterator;
    for(tmp_iterator = path.poss_perms.begin(); tmp_iterator != path.poss_perms.end(); tmp_iterator++) {
      std::vector<uint32_t> tmpPerm = (*tmp_iterator).first;
      for (uint32_t j=0; j<tmpPerm.size(); ++j)
        tmp_perms(tmpPerm[j]-1,(*tmp_iterator).second)++;
    }
    out.CreateGroup(prefix+"sectors");
    std::string data_type = "avg_pairs";
    out.Write(prefix+"sectors/data_type",data_type);
    vec<uint32_t> tmp_perm_indices(path.poss_perms.size());
    for (uint32_t i=0; i<path.poss_perms.size(); ++i)
      tmp_perm_indices(i) = i;
    out.Write(prefix+"sectors/x", tmp_perm_indices);
    out.Write(prefix+"sectors/possPerms", tmp_perms);
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
  sector_energies.clear();
}

void Energy::Accumulate()
{
  // Measure energy
  path.SetMode(1);
  double total_energy = 0.;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    if (!action_list[i]->is_importance_weight) {
      double action_energy = path.sign*path.importance_weight*action_list[i]->DActionDBeta();
      energies(i) += action_energy;
      total_energy += action_energy;
      if (measure_potential)
        potentials(i) += path.sign*path.importance_weight*action_list[i]->Potential();
    }
  }

  // Store sector with total energy
  if (measure_per_sector) {
    std::vector<uint32_t> cycle;
    path.SetCycleCount(species_i, cycle);
    uint32_t sector = path.GetPermSector(species_i, cycle);
    std::pair<uint32_t,double> sector_energy(sector,path.sign*total_energy/path.n_bead);
    sector_energies.push_back(sector_energy);
  }

  n_measure += 1;
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
      // Map out the sectors std::vector
      std::map<uint32_t,std::vector<double>> sector_map;
      for (uint32_t i=0; i<n_measure; i++) {
        std::pair<uint32_t,double> sector_energy = sector_energies.back();
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
        sector_energies.pop_back();
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
          out.CreateExtendableDataSet("/"+prefix+"/sectors/", "y", sector_info_vec);
        } else
          out.AppendDataSet("/"+prefix+"/sectors/", "y", sector_info_vec);
      }
    }

    if (first_time)
      first_time = 0;

    Reset();
  }
}
