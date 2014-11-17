#ifndef SpeciesClass_H
#define SpeciesClass_H

#include "Utils/config.h"       // Standard libraries
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"

class Species
{
private:

protected:

public:
  // Constructor
  Species(Input &in, IOClass &out, int tmpIS)
    : iS(tmpIS)
  {
    Init(in, out);
  };
  void Init(Input &in, IOClass &out);

  // Details
  int iS; // index
  string name;

  // Given Global Constants
  unsigned int nPart;
  RealType lambda;

  // Fermions
  bool fermi;
  bool fixedNode;

  // Rho_k
  bool needUpdateRhoK;
};

#endif
