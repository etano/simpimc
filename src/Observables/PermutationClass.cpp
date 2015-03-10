#include "PermutationClass.h"

void Permutation::Init(Input &in)
{
  // Read in info
  uint sectorMax = in.getAttribute<uint>("sectorMax",0);
  species = in.getAttribute<string>("species");
  out.Write("/Observables/"+name+"/species", species);
  out.Write("/Observables/"+name+"/sectorMax", sectorMax);

  // Set up permutation sectors
  path.GetSpeciesInfo(species, iSpecies);
  cycles.set_size(path.speciesList[iSpecies]->nPart);
  path.SetupPermSectors(path.speciesList[iSpecies]->nPart, sectorMax);
  firstSector = true;

  // Write out possible sectors
  mat<uint> tmpPerms;
  tmpPerms.zeros(path.speciesList[iSpecies]->nPart,path.possPerms.size());
  map<vector<uint>,uint>::iterator tmpIt;
  for(tmpIt = path.possPerms.begin(); tmpIt != path.possPerms.end(); tmpIt++) {
    vector<uint> tmpPerm = (*tmpIt).first;
    for (uint j=0; j<tmpPerm.size(); ++j)
      tmpPerms(tmpPerm[j]-1,(*tmpIt).second)++;
  }
  out.CreateGroup(prefix+"sectors");
  string data_type = "pairs";
  out.Write(prefix+"sectors/data_type",data_type);
  vec<uint> tmpPermIndices(path.possPerms.size());
  for (uint i=0; i<path.possPerms.size(); ++i)
    tmpPermIndices(i) = i;
  out.Write(prefix+"sectors/x", tmpPermIndices);
  out.Write(prefix+"sectors/possPerms", tmpPerms);

  // Write out possible cycles
  vec<uint> tmpCycles(path.speciesList[iSpecies]->nPart);
  for (uint iP=0; iP<path.speciesList[iSpecies]->nPart; ++iP)
    tmpCycles(iP) = iP+1;
  out.CreateGroup(prefix+"cycles");
  data_type = "histogram";
  out.Write(prefix+"cycles/data_type",data_type);
  out.Write(prefix+"cycles/x", tmpCycles);

  Reset();
}

void Permutation::Reset()
{
  nMeasure = 0;
  cycles.zeros();
  sectors.clear();
}

void Permutation::Accumulate()
{
  path.SetMode(1);
  vector<uint> cycle;
  path.SetCycleCount(iSpecies, cycle);
  uint sector = path.GetPermSector(iSpecies, cycle);
  sectors.push_back(sector);
  for (auto& c: cycle)
    cycles(c-1) += 1.;
  nMeasure += 1;
}

void Permutation::Write()
{
  if (nMeasure > 0) {

    // Map out the sectors vector
    map<uint,uint> sectorMap;
    for (uint i=0; i<nMeasure; i++) {
      uint sector = sectors.back();
      if (sectorMap.find(sector) == sectorMap.end())
        sectorMap.insert(pair<uint,uint>(sector,1));
      else
        sectorMap[sector]++;
      sectors.pop_back();
    }

    // Put the map into an array and write
    map<uint,uint>::iterator it;
    for(it = sectorMap.begin(); it != sectorMap.end(); it++) {
      vec<uint> tmpSectorCount(2);
      tmpSectorCount(0) = (*it).first;
      tmpSectorCount(1) = (*it).second;
      if (firstTime && firstSector) {
        firstSector = false;
        out.CreateExtendableDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
      } else
        out.AppendDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
    }

    // CycleCount
    uint nCycles = 0;
    for (uint i=0; i<cycles.size(); i++)
      nCycles += cycles(i);
    double norm = 1./nCycles;
    cycles *= norm;
    if (firstTime)
      out.CreateExtendableDataSet("/"+prefix+"/cycles/", "y", cycles);
    else
      out.AppendDataSet("/"+prefix+"/cycles/", "y", cycles);

    if (firstTime)
      firstTime = false;

    Reset();
  }
}
