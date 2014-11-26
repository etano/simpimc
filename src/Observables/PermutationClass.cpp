#include "PermutationClass.h"

void Permutation::Init(Input &in)
{
  // Read in info
  int sectorMax = in.getAttribute<int>("sectorMax",0);
  species = in.getAttribute<string>("species");
  out.Write("/Observables/"+name+"/species", species);
  out.Write("/Observables/"+name+"/sectorMax", sectorMax);

  // Set up permutation sectors
  path.GetSpeciesInfo(species, iSpecies);
  cycles.set_size(path.speciesList[iSpecies]->nPart);
  path.SetupPermSectors(path.speciesList[iSpecies]->nPart, sectorMax);
  firstSector = true;

  // Write out possible sectors
  mat<int> tmpPerms;
  tmpPerms.zeros(path.speciesList[iSpecies]->nPart,path.possPerms.size());
  map<vector<int>,int>::iterator tmpIt;
  for(tmpIt = path.possPerms.begin(); tmpIt != path.possPerms.end(); tmpIt++) {
    vector<int> tmpPerm = (*tmpIt).first;
    for (int j=0; j<tmpPerm.size(); ++j)
      tmpPerms(tmpPerm[j]-1,(*tmpIt).second)++;
  }
  out.CreateGroup(prefix+"sectors");
  string data_type = "pairs";
  out.Write(prefix+"sectors/data_type",data_type);
  vec<int> tmpPermIndices(path.possPerms.size());
  for (int i=0; i<path.possPerms.size(); ++i)
    tmpPermIndices(i) = i;
  out.Write(prefix+"sectors/x", tmpPermIndices);
  out.Write(prefix+"sectors/possPerms", tmpPerms);

  // Write out possible cycles
  vec<int> tmpCycles(path.speciesList[iSpecies]->nPart);
  for (int iP=0; iP<path.speciesList[iSpecies]->nPart; ++iP)
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
  vector<int> cycle;
  path.SetCycleCount(iSpecies, cycle);
  int sector = path.GetPermSector(iSpecies, cycle);
  sectors.push_back(sector);
  for (vector<int>::size_type j=0; j != cycle.size(); j++)
    cycles(cycle[j]-1) += 1.;
  nMeasure += 1;
}

void Permutation::Write()
{
  if (nMeasure > 0) {

    // Map out the sectors vector
    map<int,int> sectorMap;
    for (int i=0; i<nMeasure; i++) {
      int sector = sectors.back();
      if (sectorMap.find(sector) == sectorMap.end())
        sectorMap.insert(pair<int,int>(sector,1));
      else
        sectorMap[sector]++;
      sectors.pop_back();
    }

    // Put the map into an array and write
    map<int,int>::iterator it;
    for(it = sectorMap.begin(); it != sectorMap.end(); it++) {
      vec<int> tmpSectorCount(2);
      tmpSectorCount(0) = (*it).first;
      tmpSectorCount(1) = (*it).second;
      if (firstTime && firstSector) {
        firstSector = false;
        out.CreateExtendableDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
      } else
        out.AppendDataSet("/"+prefix+"/sectors/", "y", tmpSectorCount);
    }

    // CycleCount
    int nCycles = 0;
    for (int i=0; i<cycles.size(); i++)
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
