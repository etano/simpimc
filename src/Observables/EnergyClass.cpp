#include "EnergyClass.h"

void Energy::Init(Input &in)
{
  measureV = in.getAttribute<int>("measureV",0);
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
}

void Energy::Accumulate()
{
  path.SetMode(1);
  for (int i=0; i<actionList.size(); ++i) {
    if (!actionList[i]->isImportanceWeight) {
      Es(i) += path.importance_weight*actionList[i]->DActionDBeta();
      if (measureV)
        Vs(i) += path.importance_weight*actionList[i]->Potential();
    }
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

    if (firstTime)
      firstTime = 0;

    Reset();
  }
}
