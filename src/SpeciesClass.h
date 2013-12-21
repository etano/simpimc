#ifndef SpeciesClass_H
#define SpeciesClass_H

#include "config.h"       // Standard libraries
#include "IO/InputClass.h"

class Species
{
private:

protected:

public:
  // Constructor
  Species(Input &in, int nBead, int nD)
  {
    Init(in, nBead, nD);
  };
  void Init(Input &in, int nBead, int nD);

  // Details
  string name, type;

  // Given Global Constants
  unsigned int nPart;
  RealType lambda;

  // Fermions
  bool fermi;
};

#endif
