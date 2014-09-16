#ifndef BeadClass_H
#define BeadClass_H

#include "Utils/config.h"
#include "SpeciesClass.h"

class Bead
{
private:
  unsigned int nD;
protected:

public:
  // Constructor
  Bead(unsigned int tmpND, Species &tmpSpecies, unsigned int tmpP, unsigned int tmpB)
    : nD(tmpND), species(tmpSpecies), p(tmpP), b(tmpB)
  {
    self = this;
    r.zeros(nD);
    r(0) = 0.5*(p - 1.);
    storeR();
    s = species.iS;
  }

  Species &species;

  void store();
  void restore();
  void storeR();
  void restoreR();
  void storePartRecord();
  void restorePartRecord();
  void storeNodeDistance();
  void restoreNodeDistance();
  void move( Tvector& dr );
  Bead* nextB( unsigned int n );
  Bead* prevB( unsigned int n );

  unsigned int p;
  unsigned int b;
  unsigned int s;
  RealType nDist, nDistC;
  Tvector r, rC;
  Bead *self;
  Bead *next, *nextC;
  Bead *prev, *prevC;
};

#endif
