#include "PathClass.h"

// Get Kinetic Action
double Path::getK()
{
  double tot = 0.0;  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1)  {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)  {
      dr = bead(iPart,iBead) -> r - (bead(iPart,iBead) -> next -> r);
      PutInBox(dr);
      tot += dot( dr , dr );
    }
  }
  
  return oneOver4LamTau * tot;
}

// Get Single Particle Kinetic Action
double Path::getK( const int iPart )
{
  double tot = 0.0;
  dr = bead(iPart,0) -> r - (bead(iPart,0) -> prev -> r);
  PutInBox(dr);
  tot += dot( dr , dr );
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    dr = bead(iPart,iBead) -> r - (bead(iPart,iBead) -> next -> r);
    PutInBox(dr);
    tot += dot( dr , dr );
  }

  return oneOver4LamTau * tot;
}

// Get Single Bead Kinetic Action
double Path::getK( const int iPart , const int iBead )
{
  dr = bead(iPart,iBead) -> r - (bead(iPart,iBead) -> prev -> r);
  PutInBox(dr);
  double dPrev = dot( dr , dr );
  dr = bead(iPart,iBead) -> r - (bead(iPart,iBead) -> next -> r);
  PutInBox(dr);
  double dNext = dot( dr , dr );

  return oneOver4LamTau * (dPrev + dNext);
}

// Get Single Bead Kinetic Action
double Path::getK( Bead *bi )
{
  dr = bi -> r - (bi -> prev -> r);
  PutInBox(dr);
  double dPrev = dot( dr , dr );  
  dr = bi -> r - (bi -> next -> r);
  PutInBox(dr);
  double dNext = dot( dr , dr );

  return oneOver4LamTau * (dPrev + dNext);
}

// Get Potential Action
double Path::getV()
{
  double Vext = 0.0;
  double Vint = 0.0;
  for (unsigned int iPart = 0; iPart < nPart-1; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) {
        Vint += getVint( bead(iPart,iBead) , bead(jPart,iBead) );
      }
      Vext += getVext( bead(iPart,iBead) );
    }
  }
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    Vext += getVext( bead(nPart-1,iBead) );
  }

  return Vext + Vint;
}

// Get Single Particle Potential Action
double Path::getV( const unsigned int iPart )
{
  double Vext = 0.0;
  double Vint = 0.0;
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {  
    for (unsigned int jPart = 0; jPart < iPart; jPart += 1)
      Vint += getVint( bead(iPart,iBead) , bead(jPart,iBead) );
    for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) 
      Vint += getVint( bead(iPart,iBead) , bead(jPart,iBead) );

    Vext += getVext( bead(iPart,iBead) );
  }

  return  Vext + Vint;
}

// Get Single Bead Potential Action
double Path::getV( const int unsigned iPart , const int iBead )
{
  double Vext = 0.0;
  double Vint = 0.0;
  for (unsigned int jPart = 0; jPart < iPart; jPart += 1)
    Vint += getVint( bead(iPart,iBead) , bead(jPart,iBead) );
  for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) 
    Vint += getVint( bead(iPart,iBead) , bead(jPart,iBead) );

  Vext += getVext( bead(iPart,iBead) );

  return Vext + Vint;
}

// Get Single Bead Potential Action
double Path::getV( Bead *bi )
{
  double Vext = 0.0;
  double Vint = 0.0;
  const unsigned int iPart = bi -> p;
  for (unsigned int jPart = 0; jPart < iPart; jPart += 1)
    Vint += getVint( bi , bead(jPart,bi -> b) );
  for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) 
    Vint += getVint( bi , bead(jPart,bi -> b) );

  Vext += getVext( bi );

  return Vext + Vint;
}

// Get Two Bead Interaction Action
double Path::getVint( Bead *b1 , Bead *b2 )
{
  return 0;
}

// Get External Potential Action (w/ LB Correction)
double Path::getVext( Bead *b )
{
  if(trap) {
    if(mode) {
      dr = b -> r;
    } else {
      dr = b -> rC;
    }
    return halfTauOmega2 * onePlusTau2Omega2Over12 * dot( dr , dr );
  } else {
    return 0;
  }
}

// Get Single Bead Nodal Action
double Path::getN( const int iPart , const int iBead )
{
  if (!useNodeDist) return 0;
  updateNodeDistance(bead(iPart,iBead));
  updateNodeDistance(bead(iPart,iBead)->next);
  double nD1 = bead(iPart,iBead) -> nDist;
  double nD2 = bead(iPart,iBead) -> next -> nDist;
  double nD1nD2 = nD1 * nD2;
  //if (!nD1) {
  //  nD1nD2 = nD2 * nD2;
  //  cerr << "nD1" << endl;
  //} else if (!nD2) {
  //  nD1nD2 = nD1 * nD1;
  //  cerr << "nD2" << endl;
  //}
  //double factor = -log1p(-exp(-0.5*nD1nD2*oneOverLamTau));
  //double factor = -log1p(-nD1nD2*oneOverLamTau);
  double factor = -log1p(-exp(-0.5*nD1nD2*oneOverLamTau));
  return factor;
}

// Get Time Slice Nodal Action
double Path::getNSlice( const int iBead , const int skip )
{
  if (!useNodeDist) return 0;
  double N = 0;
  for (unsigned int iPart = 0; iPart < nPart; iPart++) {
    N += getN(iPart,iBead,skip);
  }
  return N;
}

// Get Single Bead Nodal Action
double Path::getN( const int iPart , const int iBead , const int skip )
{
  if (!useNodeDist) return 0;
  Bead *b1, *b2;
  b1 = bead(iPart,iBead);
  b2 = bead(iPart,iBead)->nextB(skip);
  double nD1, nD2;
  if (!mode) {
    updateNodeDistance(b1);
    updateNodeDistance(b2);
    nD1 = b1->nDistC;
    nD2 = b2->nDistC;
  } else {
    nD1 = b1->nDist;
    nD2 = b2->nDist;
  }
  if (nD1 < 0.0 || nD2 < 0.0) {
    return 1e100;
  }
  double nD1nD2 = nD1 * nD2;
  if (!nD1) {
    nD1nD2 = nD2 * nD2;
  } else if (!nD2) {
    nD1nD2 = nD1 * nD1;
  }
  double N = 0.0;
  //N += -log1p(-exp(-nD1nD2*oneOverLamTau/skip));
  N += -log1p(-exp(-0.5*nD1nD2*oneOverLamTau/skip));

  return N;
}

// Get Single Bead Nodal Action
double Path::getN( Bead *b , int skip )
{
  return getN( b->p , b->b , skip );
}

// Get Single Bead Nodal Action
double Path::getN( Bead *b )
{
  return getN( b -> p , b -> b );
}

// Get Nodal Action of Bead Set
double Path::getN( std::vector<Bead*>& beads )
{
  double N(0.0);
  for (std::vector<Bead*>::const_iterator b = beads.begin(); b != beads.end(); ++b)
    N += getN((*b));
  return N;
}
