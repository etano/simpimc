#include "EnergyClass.h"

void Energy::Init(Input &in)
{
  measureV = in.getAttribute<int>("measureV",0);

  // Initialize permutation sector things if tracking them
  measurePerSector = in.getAttribute<int>("measurePerSector",0);
  if (measurePerSector) {
    int sectorMax = in.getAttribute<int>("sectorMax",0);
    string species = in.getAttribute<string>("species");

    // Set up permutation sectors
    path.GetSpeciesInfo(species, iSpecies);
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
    string data_type = "avg_pairs";
    out.Write(prefix+"sectors/data_type",data_type);
    vec<int> tmpPermIndices(path.possPerms.size());
    for (int i=0; i<path.possPerms.size(); ++i)
      tmpPermIndices(i) = i;
    out.Write(prefix+"sectors/x", tmpPermIndices);
    out.Write(prefix+"sectors/possPerms", tmpPerms);
  }

  // Resize things
  Es.set_size(actionList.size());
  if (measureV)
    Vs.set_size(actionList.size());
  Reset();
}

void Energy::Reset()
{
  nMeasure = 0;
  Es.zeros();
  if (measureV)
    Vs.zeros();
  sectorEs.clear();
}

void Energy::Accumulate()
{
  // Measure energy
  path.SetMode(1);
  double totE = 0.;
  for (int i=0; i<actionList.size(); ++i) {
    if (!actionList[i]->isImportanceWeight) {
      double actionE = path.sign*path.importance_weight*actionList[i]->DActionDBeta();
      Es(i) += actionE;
      totE += actionE;
      if (measureV)
        Vs(i) += path.sign*path.importance_weight*actionList[i]->Potential();
    }
  }

  // Store sector with total energy
  if (measurePerSector) {
    vector<int> cycle;
    path.SetCycleCount(iSpecies, cycle);
    int sector = path.GetPermSector(iSpecies, cycle);
    std::pair<int,double> sectorE(sector,path.sign*totE/path.nBead);
    sectorEs.push_back(sectorE);
  }

  nMeasure += 1;
}

void Energy::Write()
{
  if (nMeasure > 0) {
    double norm = path.nBead*nMeasure;

    // Write Es
    Es = Es/norm;
    double E = sum(Es);
    if (firstTime) {
      out.CreateGroup(prefix+"Total");
      out.CreateExtendableDataSet("/"+prefix+"Total/", "x", E);
      string data_type = "scalar";
      out.Write(prefix+"Total/data_type",data_type);
      for (int i=0; i<actionList.size(); ++i) {
        out.CreateGroup(prefix+actionList[i]->name);
        out.CreateExtendableDataSet("/"+prefix+actionList[i]->name+"/", "x", Es(i));
        out.Write(prefix+actionList[i]->name+"/data_type", data_type);
      }
    } else {
      out.AppendDataSet("/"+prefix+"Total/", "x", E);
      for (int i=0; i<actionList.size(); ++i)
        out.AppendDataSet("/"+prefix+actionList[i]->name+"/", "x", Es(i));
    }

    // Write Vs
    if (measureV) {
      Vs = Vs/norm;
      double V = sum(Vs);
      if (firstTime) {
        out.CreateGroup(prefix+"VTotal");
        out.CreateExtendableDataSet("/"+prefix+"VTotal/", "x", V);
        string data_type = "scalar";
        out.Write(prefix+"VTotal/data_type",data_type);
        for (int i=0; i<actionList.size(); ++i) {
          out.CreateGroup(prefix+"V"+actionList[i]->name);
          out.CreateExtendableDataSet("/"+prefix+"V"+actionList[i]->name+"/", "x", Vs(i));
          out.Write(prefix+"V"+actionList[i]->name+"/data_type", data_type);
        }
      } else {
        out.AppendDataSet("/"+prefix+"VTotal/", "x", V);
        for (int i=0; i<actionList.size(); ++i)
          out.AppendDataSet("/"+prefix+"V"+actionList[i]->name+"/", "x", Vs(i));
      }
    }

    // Write sector Es
    if (measurePerSector) {
      // Map out the sectors vector
      std::map<int,std::vector<double> > sectorMap;
      for (int i=0; i<nMeasure; i++) {
        std::pair<int,double> sectorE = sectorEs.back();
        int sector = sectorE.first;
        double energy = sectorE.second;
        if (sectorMap.find(sector) == sectorMap.end()) {
          std::vector<double> sectorInfo = {energy,0.,1};
          sectorMap.insert(std::pair<int,std::vector<double> >(sector, sectorInfo));
        } else {
          double E0 = sectorMap[sector][0];
          double var0 = sectorMap[sector][1];
          double N0 = sectorMap[sector][2];
          double N1 = N0 + 1;
          double E1 = (N0*E0 + energy)/N1; // New average
          double var1 = ((N0*var0 + N0*E0*E0 + energy*energy)/N1) - E1*E1;
          sectorMap[sector][0] = E1;
          sectorMap[sector][1] = var1;
          sectorMap[sector][2] = N1;
        }
        sectorEs.pop_back();
      }

      // Put the map into an array and write
      map<int,int>::iterator it;
      for(auto& sectorInfo: sectorMap) {
        vec<double> sectorInfoVec(4);
        sectorInfoVec(0) = sectorInfo.first;
        sectorInfoVec(1) = sectorInfo.second[0];
        sectorInfoVec(2) = sectorInfo.second[1];
        sectorInfoVec(3) = sectorInfo.second[2];
        if (firstTime && firstSector) {
          firstSector = false;
          out.CreateExtendableDataSet("/"+prefix+"/sectors/", "y", sectorInfoVec);
        } else
          out.AppendDataSet("/"+prefix+"/sectors/", "y", sectorInfoVec);
      }
    }

    if (firstTime)
      firstTime = 0;

    Reset();
  }
}
