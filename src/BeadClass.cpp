#include "BeadClass.h"

void Bead::store()
{
  storeR();
  storeRhoK();
  storePartRecord();
  storeNodeDistance();
}

void Bead::restore()
{
  restoreR();
  restoreRhoK();
  restorePartRecord();
  restoreNodeDistance();
}

void Bead::storePartRecord()
{
  nextC = next;
  prevC = prev;
}

void Bead::restorePartRecord()
{
  next = nextC;
  prev = prevC;
}

Bead* Bead::nextB( unsigned int n )
{
  Bead *bead = this;
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> next;
  return bead;
}

Bead* Bead::nextBC( unsigned int n )
{
  Bead *bead = this;
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> nextC;
  return bead;
}

Bead* Bead::prevB( unsigned int n )
{
  Bead *bead = this;
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> prev;
  return bead;
}

Bead* Bead::prevBC( unsigned int n )
{
  Bead *bead = this;
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> prevC;
  return bead;
}

