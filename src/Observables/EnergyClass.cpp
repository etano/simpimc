#include "EnergyClass.h"

void Energy::Init()
{
  Es.resize(actionList.size());
  Reset();
}

void Energy::Reset()
{
  nMeasure = 0;
  for (int i=0; i<Es.size(); ++i)
    Es[i] = 0.;
}

void Energy::Accumulate()
{
  for (int i=0; i<actionList.size(); ++i)
    Es[i] += actionList[i]->DActionDBeta();
  nMeasure += 1;
}

void Energy::Write()
{
  RealType norm = 1./(path.nBead*nMeasure);
  RealType E = 0.;
  for (int i=0; i<actionList.size(); ++i) {
    Es[i] *= norm;
    E += Es[i];
  }

  if (FirstTime) {
    FirstTime = 0;
    out.CreateExtendableDataSet("/"+name+"/", "Total", E);
    for (int i=0; i<actionList.size(); ++i)
      out.CreateExtendableDataSet("/"+name+"/", actionList[i]->name, Es[i]);
  } else {
    out.AppendDataSet("/"+name+"/", "Total", E);
    for (int i=0; i<actionList.size(); ++i)
      out.AppendDataSet("/"+name+"/", actionList[i]->name, Es[i]);
  }

  Reset();
}
