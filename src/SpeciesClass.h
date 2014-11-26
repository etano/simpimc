#ifndef SpeciesClass_H
#define SpeciesClass_H

#include "config.h"
#include "BeadClass.h"

class Species
{
private:

protected:

public:
  // Constructor
  Species(Input &in, IOClass &out, int t_iS, int t_nD, int t_nBead)
    : iS(t_iS), nD(t_nD), nBead(t_nBead)
  {
    Init(in, out);
  };
  void Init(Input &in, IOClass &out);
  void InitPaths(Input &in, IOClass &out, RNG &rng, CommunicatorClass& InterComm, int L);

  // Beads
  field<Bead*> bead;
  string initType;

  // Given Global Constants
  unsigned int nPart, nD, nBead, iS;
  string name;
  double lambda;

  // Fermions
  bool fermi;
  bool fixedNode;

  // Rho_k
  bool needUpdateRhoK;
};

#endif
