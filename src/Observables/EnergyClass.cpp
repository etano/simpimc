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
    Es(i) += actionList[i]->DActionDBeta();
    if (measureV)
      Vs(i) += actionList[i]->Potential();
  }
  nMeasure += 1;
}

void Energy::Write()
{
  if (nMeasure > 0) {
    RealType norm = path.nBead*nMeasure;

    // Write Es
    Es = Es/norm;
    RealType E = sum(Es);
    if (firstTime) {
      out.CreateExtendableDataSet("/Observables/"+name+"/", "Total", E);
      for (int i=0; i<actionList.size(); ++i)
        out.CreateExtendableDataSet("/Observables/"+name+"/", actionList[i]->name, Es(i));
    } else {
      out.AppendDataSet("/Observables/"+name+"/", "Total", E);
      for (int i=0; i<actionList.size(); ++i)
        out.AppendDataSet("/Observables/"+name+"/", actionList[i]->name, Es(i));
    }

    // Write Vs
    if (measureV) {
      Vs = Vs/norm;
      RealType V = sum(Vs);
      if (firstTime) {
        out.CreateExtendableDataSet("/Observables/"+name+"/", "vTotal", V);
        for (int i=0; i<actionList.size(); ++i)
          out.CreateExtendableDataSet("/Observables/"+name+"/", "v"+actionList[i]->name, Vs(i));
      } else {
        out.AppendDataSet("/Observables/"+name+"/", "vTotal", V);
        for (int i=0; i<actionList.size(); ++i)
         out.AppendDataSet("/Observables/"+name+"/", "v"+actionList[i]->name, Vs(i));
      }
    }

    if (firstTime)
      firstTime = 0;

    Reset();
  }
}
