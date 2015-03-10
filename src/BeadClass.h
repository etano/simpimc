#ifndef BeadClass_H
#define BeadClass_H

#include "config.h"

struct Bead
{
public:
  Bead();
  Bead(uint tmpND, uint tmpS, uint tmpP, uint tmpB)
    : nD(tmpND), s(tmpS), p(tmpP), b(tmpB), r(tmpND), rC(tmpND), isIra(false), isMasha(false)
  {}

  uint p, b, s, nD;
  bool isIra, isMasha; // head and tail (respectively)
  double nDist, nDistC;
  vec<double> r, rC;
  vec<complex<double> > rhoK, rhoKC;
  std::shared_ptr<Bead> next, nextC, prev, prevC;

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
  inline void move(const vec<double>& dr) { r += dr; };

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

  std::shared_ptr<Bead> nextB(const uint n)
  {
    std::shared_ptr<Bead> bead(next);
    for (uint i=1; i<n; i++) bead = bead -> next;
    return bead;
  }

  std::shared_ptr<Bead> nextBC(const uint n)
  {
    std::shared_ptr<Bead> bead(nextC);
    for (uint i=1; i<n; i++) bead = bead -> nextC;
    return bead;
  }

  std::shared_ptr<Bead> prevB(const uint n)
  {
    std::shared_ptr<Bead> bead(prev);
    for (uint i=1; i<n; i++) bead = bead -> prev;
    return bead;
  }

  std::shared_ptr<Bead> prevBC(const uint n)
  {
    std::shared_ptr<Bead> bead(prevC);
    for (uint i=1; i<n; i++) bead = bead -> prevC;
    return bead;
  }

};

#endif
