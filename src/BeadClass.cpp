#include "BeadClass.h"

void Bead::store()
{
  storeR();
  storePartRecord();
  storeNodeDistance();
}

void Bead::restore()
{
  restoreR();
  restorePartRecord();
  restoreNodeDistance();
}

void Bead::storeR()
{
  rC = r;
}

void Bead::restoreR()
{
  r = rC;
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

void Bead::storeNodeDistance()
{
  nDistC = nDist;
}

void Bead::restoreNodeDistance()
{
  nDist = nDistC;
}

void Bead::move( vec& dr )
{
  r += dr;
}

Bead* Bead::nextB( unsigned int n )
{
  Bead *bead = this;
  
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> next;
  
  return bead;
}

Bead* Bead::prevB( unsigned int n )
{
  Bead *bead = this;
  
  for (unsigned int i = 0; i < n; i += 1) bead = bead -> prev;
  
  return bead;
}

