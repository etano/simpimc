#include "SpeciesClass.h"

void Species::Init(Input &in, int nBead, int nD)
{
  name = in.getAttribute<string>("name");
  type = in.getAttribute<string>("type");
  nPart = in.getAttribute<int>("nPart");
  lambda = in.getAttribute<RealType>("lambda");
  fermi = in.getAttribute<int>("fermi", 0);
}
