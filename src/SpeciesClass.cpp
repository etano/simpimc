#include "SpeciesClass.h"

void Species::Init(Input &in, int nBead, int nD)
{
  name = in.get<string>("Name");
  type = in.get<string>("Type");
  nPart = in.get<int>("nPart");
  lambda = in.get<RealType>("lambda");
  fermi = in.get<int>("fermi", 0);
}
