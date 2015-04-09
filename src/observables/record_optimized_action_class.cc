#include "record_optimized_action_class.h"

void RecordOptimizedAction::Init(Input &in)
{
  // Read in action name
  std::string action_name = in.GetAttribute<std::string>("action_name");

  // Select action from list
  for (auto& t_action : action_list) {
    if (t_action->name == action_name)
      action = t_action;
  }

  // Set up param_set_count
  param_set_count.set_size(action->GetNumParamSets());
  vec<uint32_t> param_set_indices(param_set_count.size());
  for (uint32_t i=0; i<param_set_count.size(); ++i)
    param_set_indices(i) = i;
  out.Write(prefix+"x",param_set_indices);

  Reset();
}

void RecordOptimizedAction::Reset()
{
  n_measure = 0;
  param_set_count.zeros();
}

void RecordOptimizedAction::Accumulate()
{
  path.SetMode(NEW_MODE);
  param_set_count(action->GetParamSet()) += 1.;
  n_measure++;
}

void RecordOptimizedAction::Write()
{
  if (n_measure > 0) {
    double norm = 1.*n_measure;

    // Write sign
    param_set_count /= norm;
    if (first_time) {
      first_time = 0;
      out.CreateExtendableDataSet(prefix, "y", param_set_count);
    } else {
      out.AppendDataSet(prefix, "y", param_set_count);
    }

    Reset();
  }
}
