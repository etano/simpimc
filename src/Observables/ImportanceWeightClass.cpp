#include "ImportanceWeightClass.h"

void ImportanceWeight::Init(Input &in)
{
  IWs.set_size(actionList.size());
  Reset();
}

void ImportanceWeight::Reset()
{
  nMeasure = 0;
  IWs.zeros();
}

void ImportanceWeight::Accumulate()
{
  path.SetMode(1);
  path.importance_weight = path.sign;
  for (int i=0; i<actionList.size(); ++i) {
    double tmpIW = 1.;
    if (actionList[i]->isImportanceWeight)
      tmpIW = actionList[i]->ImportanceWeight();
    IWs(i) += tmpIW;
    path.importance_weight *= tmpIW;
  }
  nMeasure += 1;
}

void ImportanceWeight::Write()
{
  if (nMeasure > 0) {
    double norm = nMeasure;

    // Write IWs
    IWs = IWs/norm;
    double IW = prod(IWs);
    if (firstTime) {
      out.CreateGroup(prefix+"Total");
      out.CreateExtendableDataSet("/"+prefix+"Total/", "x", IW);
      string data_type = "scalar";
      out.Write(prefix+"Total/data_type",data_type);
      for (int i=0; i<actionList.size(); ++i) {
        out.CreateGroup(prefix+actionList[i]->name);
        out.CreateExtendableDataSet("/"+prefix+actionList[i]->name+"/", "x", IWs(i));
        out.Write(prefix+actionList[i]->name+"/data_type", data_type);
      }
      firstTime = 0;
    } else {
      out.AppendDataSet("/"+prefix+"Total/", "x", IW);
      for (int i=0; i<actionList.size(); ++i)
        out.AppendDataSet("/"+prefix+actionList[i]->name+"/", "x", IWs(i));
    }

    Reset();
  }
}
