#include "SignClass.h"

void Sign::Init(Input &in)
{
  string data_type = "scalar";
  out.Write(prefix+"/data_type",data_type);

  Reset();
}

void Sign::Reset()
{
  nMeasure = 0;
  sign = 0.;
}

void Sign::Accumulate()
{
  path.SetMode(1);
  sign += path.sign;
  nMeasure += 1;
}

void Sign::Write()
{
  if (nMeasure > 0) {
    RealType norm = 1.*nMeasure;

    // Write sign
    sign /= norm;
    if (firstTime) {
      firstTime = 0;
      out.CreateExtendableDataSet("/Observables/"+name+"/", "x", sign);
    } else {
      out.AppendDataSet("/Observables/"+name+"/", "x", sign);
    }

    Reset();
  }
}
