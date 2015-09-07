#ifndef SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PERMUTATION_CLASS_H_

#include "observable_class.h"

/// Tracks the permutation sectors and cycles
class Permutation : public Observable
{
private:
  bool first_sector; ///< Whether or not the first permutation sector has been written
  std::shared_ptr<Species> species; ///< Pointer to the relevant species
  std::vector<uint32_t> sectors; ///< Vector of sector indices
  vec<double> cycles; ///< Vector of cycle counts

  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    auto cycle = species->GetCycleCount();
    sectors.push_back(species->GetPermSector(cycle));
    for (auto& c: cycle)
      cycles(c-1) += 1;
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
    std::string species_name = in.GetAttribute<std::string>("species");
    out.Write("/Observables/"+name+"/species", species_name);
    out.Write("/Observables/"+name+"/sector_max", sector_max);

    // Set up permutation sectors
    species = path.GetSpecies(species_name);
    cycles.set_size(species->GetNPart());
    species->SetupPermSectors(sector_max);
    first_sector = true;

    // Write out possible sectors
    mat<uint32_t> perms(zeros<mat<uint32_t>>(species->GetNPart(),species->GetNPerm()));
    for(const auto& poss_perm : species->GetPossPerms()) {
      std::vector<uint32_t> perm = poss_perm.first;
      for (uint32_t j=0; j<perm.size(); ++j)
        perms(perm[j]-1,poss_perm.second)++;
    }
    out.CreateGroup(prefix+"sectors");
    std::string data_type = "pairs";
    out.Write(prefix+"sectors/data_type",data_type);
    vec<uint32_t> perm_indices(species->GetNPerm());
    for (uint32_t i=0; i<species->GetNPerm(); ++i)
      perm_indices(i) = i;
    out.Write(prefix+"sectors/x", perm_indices);
    out.Write(prefix+"sectors/possPerms", perms);

    // Write out possible cycles
    vec<uint32_t> t_cycles(species->GetNPart());
    for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
      t_cycles(p_i) = p_i+1;
    out.CreateGroup(prefix+"cycles");
    data_type = "histogram";
    out.Write(prefix+"cycles/data_type",data_type);
    out.Write(prefix+"cycles/x", t_cycles);

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
