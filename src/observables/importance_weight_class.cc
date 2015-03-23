#include "importance_weight_class.h"

void ImportanceWeight::Init(Input &in)
{
  IWs.set_size(action_list.size());
  Reset();
}

void ImportanceWeight::Reset()
{
  n_measure = 0;
  IWs.zeros();
}

void ImportanceWeight::Accumulate()
{
  path.SetMode(NEW_MODE);
  path.importance_weight = path.sign;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    double tmpIW = 1.;
    if (action_list[i]->is_importance_weight)
      tmpIW = action_list[i]->ImportanceWeight();
    IWs(i) += tmpIW;
    path.importance_weight *= tmpIW;
  }
  n_measure += 1;
}

void ImportanceWeight::Write()
{
  if (n_measure > 0) {
    double norm = n_measure;

    // Write IWs
    IWs = IWs/norm;
    double IW = prod(IWs);
    if (first_time) {
      out.CreateGroup(prefix+"Total");
      out.CreateExtendableDataSet("/"+prefix+"Total/", "x", IW);
      std::string data_type = "scalar";
      out.Write(prefix+"Total/data_type",data_type);
      for (uint32_t i=0; i<action_list.size(); ++i) {
        out.CreateGroup(prefix+action_list[i]->name);
        out.CreateExtendableDataSet("/"+prefix+action_list[i]->name+"/", "x", IWs(i));
        out.Write(prefix+action_list[i]->name+"/data_type", data_type);
      }
      first_time = 0;
    } else {
      out.AppendDataSet("/"+prefix+"Total/", "x", IW);
      for (uint32_t i=0; i<action_list.size(); ++i)
        out.AppendDataSet("/"+prefix+action_list[i]->name+"/", "x", IWs(i));
    }

    Reset();
  }
}
