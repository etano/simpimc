#ifndef SpeciesClass_H
#define SpeciesClass_H

#include "Utils/config.h"       // Standard libraries
#include "Utils/IO/InputClass.h"

class Species
{
private:

protected:

public:
  // Constructor
  Species(Input &in, int tmpIS, int nBead, int nD)
    : iS(tmpIS)
  {
    Init(in, nBead, nD);
  };
  void Init(Input &in, int nBead, int nD);

  // Details
  int iS; // index
  string name, type;

  // Given Global Constants
  unsigned int nPart;
  RealType lambda;

  // Fermions
  bool fermi;
};

#endif
