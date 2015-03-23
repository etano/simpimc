#include "sign_class.h"

void Sign::Init(Input &in)
{
  std::string data_type = "scalar";
  out.Write(prefix+"data_type",data_type);

  Reset();
}

void Sign::Reset()
{
  n_measure = 0;
  sign = 0.;
}

void Sign::Accumulate()
{
  path.SetMode(NEW_MODE);
  sign += path.sign;
  n_measure += 1;
}

void Sign::Write()
{
  if (n_measure > 0) {
    double norm = 1.*n_measure;

    // Write sign
    sign /= norm;
    if (first_time) {
      first_time = 0;
      out.CreateExtendableDataSet(prefix, "x", sign);
    } else {
      out.AppendDataSet(prefix, "x", sign);
    }

    Reset();
  }
}
