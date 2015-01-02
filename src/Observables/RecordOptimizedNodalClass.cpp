#include "RecordOptimizedNodalClass.h"

void RecordOptimizedNodal::Init(Input &in)
{
  // Read in action name
  string actionName = in.getAttribute<string>("actionName");

  // Select action from list
  for (auto& t_action : actionList) {
    if (t_action->name == actionName)
      action = t_action;
  }

  // Set up paramSetCount
  paramSetCount.set_size(action->GetNumParamSets());
  vec<int> paramSetIndices(paramSetCount.size());
  for (int i=0; i<paramSetCount.size(); ++i)
    paramSetIndices(i) = i;
  out.Write(prefix+"x",paramSetIndices);

  // Write out data type
  string data_type = "histogram";
  out.Write(prefix+"data_type",data_type);

  Reset();
}

void RecordOptimizedNodal::Reset()
{
  nMeasure = 0;
  paramSetCount.zeros();
}

void RecordOptimizedNodal::Accumulate()
{
  path.SetMode(1);
  paramSetCount(action->GetParamSet()) += 1.;
  nMeasure++;
}

void RecordOptimizedNodal::Write()
{
  if (nMeasure > 0) {
    double norm = 1.*nMeasure;

    // Write sign
    paramSetCount /= norm;
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet(prefix, "y", paramSetCount);
    } else {
      out.AppendDataSet(prefix, "y", paramSetCount);
    }

    Reset();
  }
}
