#include "SpeciesClass.h"

void Species::Init(Input &in, IOClass &out, int nBead, int nD)
{
  // Read inputs
  name = in.getAttribute<string>("name");
  type = in.getAttribute<string>("type");
  nPart = in.getAttribute<int>("nPart");
  lambda = in.getAttribute<RealType>("lambda");
  fermi = in.getAttribute<int>("fermi", 0);
  fixedNode = in.getAttribute<int>("fixedNode", 0);

  // Write to file
  out.CreateGroup("System/Particles/"+name);
  out.Write("System/Particles/"+name+"/type",type);
  out.Write("System/Particles/"+name+"/nPart",nPart);
  out.Write("System/Particles/"+name+"/lambda",lambda);
  out.Write("System/Particles/"+name+"/fermi",fermi);
  out.Write("System/Particles/"+name+"/fixedNode",fixedNode);

  // Set defaults
  needUpdateRhoK = true;
}
