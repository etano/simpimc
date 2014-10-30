#include "PermutationClass.h"

void Permutation::Init(Input &in)
{
  // Read in info
  int sectorMax = in.getAttribute<int>("sectorMax",0);
  species = in.getAttribute<string>("species");
  out.Write("/Observables/"+name+"/species", species);
  out.Write("/Observables/"+name+"/sectorMax", sectorMax);

  // Set up permutation sectors
  path.GetSpeciesInfo(species, iSpecies, offset);
  cycles.set_size(path.speciesList[iSpecies]->nPart);
  path.SetupPermSectors(path.speciesList[iSpecies]->nPart, sectorMax);
  firstSector = true;

  // Write out possible sectors
  Imatrix tmpPerms;
  tmpPerms.zeros(path.speciesList[iSpecies]->nPart,path.possPerms.size());
  map<vector<int>,int>::iterator tmpIt;
  for(tmpIt = path.possPerms.begin(); tmpIt != path.possPerms.end(); tmpIt++) {
    vector<int> tmpPerm = (*tmpIt).first;
    for (int j=0; j<tmpPerm.size(); ++j)
      tmpPerms(tmpPerm[j]-1,(*tmpIt).second)++;
  }
  out.Write("/Observables/"+name+"/poss_sectors", tmpPerms);

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
      Ivector tmpSectorCount(2);
      tmpSectorCount(0) = (*it).first;
      tmpSectorCount(1) = (*it).second;
      if (firstTime && firstSector) {
        firstSector = false;
        out.CreateExtendableDataSet("/Observables/"+name+"/", "sectors", tmpSectorCount);
      } else
        out.AppendDataSet("/Observables/"+name+"/", "sectors", tmpSectorCount);
    }

    // CycleCount
    int nCycles = 0;
    for (int i=0; i<cycles.size(); i++)
      nCycles += cycles(i);
    RealType norm = 1./nCycles;
    cycles *= norm;
    if (firstTime)
      out.CreateExtendableDataSet("/Observables/"+name+"/", "cycles", cycles);
    else
      out.AppendDataSet("/Observables/"+name+"/", "cycles", cycles);

    if (firstTime)
      firstTime = false;

    Reset();
  }
}
