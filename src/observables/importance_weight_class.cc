#include "importance_weight_class.h"

void ImportanceWeight::Init(Input &in)
{
  importance_weights.set_size(action_list.size());
  Reset();
}

void ImportanceWeight::Reset()
{
  n_measure = 0;
  importance_weights.zeros();
}

void ImportanceWeight::Accumulate()
{
  path.SetMode(NEW_MODE);
  path.importance_weight = path.sign;
  for (uint32_t i=0; i<action_list.size(); ++i) {
    double tmp_iw = 1.;
    if (action_list[i]->is_importance_weight)
      tmp_iw = action_list[i]->ImportanceWeight();
    importance_weights(i) += tmp_iw;
    path.importance_weight *= tmp_iw;
  }
  n_measure += 1;
}

void ImportanceWeight::Write()
{
  if (n_measure > 0) {
    double norm = n_measure;

    // Write importance_weights
    importance_weights = importance_weights/norm;
    double IW = prod(importance_weights);
    if (first_time) {
      out.CreateGroup(prefix+"Total");
      out.CreateExtendableDataSet("/"+prefix+"Total/", "x", IW);
      std::string data_type = "scalar";
      out.Write(prefix+"Total/data_type",data_type);
      for (uint32_t i=0; i<action_list.size(); ++i) {
        out.CreateGroup(prefix+action_list[i]->name);
        out.CreateExtendableDataSet("/"+prefix+action_list[i]->name+"/", "x", importance_weights(i));
        out.Write(prefix+action_list[i]->name+"/data_type", data_type);
      }
      first_time = 0;
    } else {
      out.AppendDataSet("/"+prefix+"Total/", "x", IW);
      for (uint32_t i=0; i<action_list.size(); ++i)
        out.AppendDataSet("/"+prefix+action_list[i]->name+"/", "x", importance_weights(i));
    }

    Reset();
  }
}
