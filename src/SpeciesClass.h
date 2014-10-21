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
  Species(Input &in, IOClass &out, int tmpIS, int nBead, int nD)
    : iS(tmpIS)
  {
    Init(in, out, nBead, nD);
  };
  void Init(Input &in, IOClass &out, int nBead, int nD);

  // Details
  int iS; // index
  string name, type;

  // Given Global Constants
  unsigned int nPart;
  RealType lambda;

  // Fermions
  bool fermi;

  // Rho_k
  bool needUpdateRhoK;
};

#endif
