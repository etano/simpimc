#ifndef BeadClass_H
#define BeadClass_H

#include "config.h"

struct Bead
{
public:
  Bead();
  ~Bead();
  Bead(unsigned int tmpND, int tmpS, unsigned int tmpP, unsigned int tmpB, unsigned int tmpId)
    : nD(tmpND), s(tmpS), p(tmpP), b(tmpB), id(tmpId), r(tmpND), rC(tmpND)
  {}

  unsigned int p, b, s, id, nD;
  double nDist, nDistC;
  vec<double> r, rC;
  vec< complex<double> > rhoK, rhoKC;
  Bead *next, *nextC, *prev, *prevC;

  inline void storeR() { rC = r; };
  inline void restoreR() { r = rC; };
  inline void storeRhoK() { rhoKC = rhoK; };
  inline void restoreRhoK() { rhoK = rhoKC; };
  inline void storePrev() { prevC = prev; };
  inline void restorePrev() { prev = prevC; };
  inline void storeNext() { nextC = next; };
  inline void restoreNext() { next = nextC; };
  inline void storeNodeDistance() { nDistC = nDist; };
  inline void restoreNodeDistance() { nDist = nDistC; };
  inline void move( vec<double>& dr ) { r += dr; };

  inline void store()
  {
    storeR();
    storeRhoK();
    storePartRecord();
    storeNodeDistance();
  }

  inline void restore()
  {
    restoreR();
    restoreRhoK();
    restorePartRecord();
    restoreNodeDistance();
  }

  inline void storePartRecord()
  {
    nextC = next;
    prevC = prev;
  }

  inline void restorePartRecord()
  {
    next = nextC;
    prev = prevC;
  }

  Bead* nextB(unsigned int const n)
  {
    Bead *bead = this;
    for (unsigned int i=0; i<n; i++) bead = bead -> next;
    return bead;
  }

  Bead* nextBC(unsigned int const n)
  {
    Bead *bead = this;
    for (unsigned int i=0; i<n; i++) bead = bead -> nextC;
    return bead;
  }

  Bead* prevB(unsigned int const n)
  {
    Bead *bead = this;
    for (unsigned int i=0; i<n; i++) bead = bead -> prev;
    return bead;
  }

  Bead* prevBC(unsigned int const n)
  {
    Bead *bead = this;
    for (unsigned int i=0; i<n; i++) bead = bead -> prevC;
    return bead;
  }

};

#endif
