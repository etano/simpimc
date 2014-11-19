#ifndef SpeciesClass_H
#define SpeciesClass_H

#include "Utils/config.h"       // Standard libraries
#include "Utils/Communication/Communication.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"
#include "BeadClass.h"

class Species
{
private:

protected:

public:
  // Constructor
  Species(Input &in, IOClass &out, int t_iS, int t_offset, int t_nD, int t_nBead)
    : iS(t_iS), offset(t_offset), nD(t_nD), nBead(t_nBead)
  {
    Init(in, out);
  };
  void Init(Input &in, IOClass &out);
  void InitPaths(Input &in, IOClass &out, RNG &rng, field<Bead*>& bead, CommunicatorClass& InterComm, int L);
  string initType;

  // Given Global Constants
  unsigned int nPart, nD, nBead, iS, offset;
  string name;
  RealType lambda;

  // Fermions
  bool fermi;
  bool fixedNode;

  // Rho_k
  bool needUpdateRhoK;
};

#endif
