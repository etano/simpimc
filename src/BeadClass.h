#ifndef BeadClass_H
#define BeadClass_H

#include "StandardLibraries.h"       // Standard libraries
#include "GlobalConstants.h"

class Bead
{
private:   
  unsigned int nD;
protected:
  
public:
  // Constructor
  Bead( unsigned int nDIn , unsigned int pIn , unsigned int bIn )
    : nD(nDIn) , p(pIn) , b(bIn)
  {
    self = this;
    r.zeros(nD);
    r.fill(0.01*p);
    storeR();
  }  
  
  void store();
  void restore();
  void storeR();
  void restoreR();
  void storePartRecord();
  void restorePartRecord();
  void storeNodeDistance();
  void restoreNodeDistance();
  void move( vec& dr );
  Bead* nextB( unsigned int n );
  Bead* prevB( unsigned int n );

  unsigned int p;
  unsigned int b;
  double nDist, nDistC;  
  vec r, rC;
  Bead *self;
  Bead *next, *nextC;
  Bead *prev, *prevC;
};

#endif
