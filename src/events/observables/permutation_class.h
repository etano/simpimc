#ifndef SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_

#include "observable_class.h"

/// Tracks the permutation sectors and cycles
class Permutation : public Observable
{
private:
  bool first_sector; ///< Whether or not the first permutation sector has been written
  uint32_t species_i; ///< Index of relevant species
  std::string species; ///< Name of relevant species
  std::vector<uint32_t> sectors; ///< Vector of sector indices
  vec<double> cycles; ///< Vector of cycle counts

  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    std::vector<uint32_t> cycle;
    path.SetCycleCount(species_i, cycle);
    uint32_t sector = path.GetPermSector(species_i, cycle);
    sectors.push_back(sector);
    for (auto& c: cycle)
      cycles(c-1) += 1.;
    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_measure = 0;
    cycles.zeros();
    sectors.clear();
  }

public:
  /// Constructor calls Init
  Permutation(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    // Read in info
    uint32_t sector_max = in.GetAttribute<uint32_t>("sector_max",0);
    species = in.GetAttribute<std::string>("species");
    out.Write("/Observables/"+name+"/species", species);
    out.Write("/Observables/"+name+"/sector_max", sector_max);

    // Set up permutation sectors
    path.GetSpeciesInfo(species, species_i);
    cycles.set_size(path.species_list[species_i]->n_part);
    path.SetupPermSectors(species_i, sector_max);
    first_sector = true;

    // Write out possible sectors
    mat<uint32_t> tmp_perms;
    tmp_perms.zeros(path.species_list[species_i]->n_part,path.poss_perms[species_i].size());
    std::map<std::vector<uint32_t>,uint32_t>::iterator tmp_iterator;
    for(tmp_iterator = path.poss_perms[species_i].begin(); tmp_iterator != path.poss_perms[species_i].end(); tmp_iterator++) {
      std::vector<uint32_t> tmpPerm = (*tmp_iterator).first;
      for (uint32_t j=0; j<tmpPerm.size(); ++j)
        tmp_perms(tmpPerm[j]-1,(*tmp_iterator).second)++;
    }
    out.CreateGroup(prefix+"sectors");
    std::string data_type = "pairs";
    out.Write(prefix+"sectors/data_type",data_type);
    vec<uint32_t> tmp_perm_indices(path.poss_perms[species_i].size());
    for (uint32_t i=0; i<path.poss_perms[species_i].size(); ++i)
      tmp_perm_indices(i) = i;
    out.Write(prefix+"sectors/x", tmp_perm_indices);
    out.Write(prefix+"sectors/possPerms", tmp_perms);

    // Write out possible cycles
    vec<uint32_t> tmp_cycles(path.species_list[species_i]->n_part);
    for (uint32_t p_i=0; p_i<path.species_list[species_i]->n_part; ++p_i)
      tmp_cycles(p_i) = p_i+1;
    out.CreateGroup(prefix+"cycles");
    data_type = "histogram";
    out.Write(prefix+"cycles/data_type",data_type);
    out.Write(prefix+"cycles/x", tmp_cycles);

    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {

      // Map out the sectors std::vector
      std::map<uint32_t,uint32_t> sectorMap;
      for (uint32_t i=0; i<n_measure; i++) {
        uint32_t sector = sectors.back();
        if (sectorMap.find(sector) == sectorMap.end())
          sectorMap.insert(std::pair<uint32_t,uint32_t>(sector,1));
        else
          sectorMap[sector]++;
        sectors.pop_back();
      }

      // Put the std::map into an array and write
      std::map<uint32_t,uint32_t>::iterator it;
      for(it = sectorMap.begin(); it != sectorMap.end(); it++) {
        vec<uint32_t> tmpSectorCount(2);
        tmpSectorCount(0) = (*it).first;
        tmpSectorCount(1) = (*it).second;
        if (first_time && first_sector) {
          first_sector = false;
          out.CreateExtendableDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
        } else
          out.AppendDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
      }

      // CycleCount
      uint32_t nCycles = 0;
      for (uint32_t i=0; i<cycles.size(); i++)
        nCycles += cycles(i);
      double norm = 1./nCycles;
      cycles *= norm;
      if (first_time)
        out.CreateExtendableDataSet("/"+prefix+"/cycles/", "y", cycles);
      else
        out.AppendDataSet("/"+prefix+"/cycles/", "y", cycles);

      if (first_time)
        first_time = false;

      Reset();
    }
  }

};

#endif // SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
