#include "PathClass.h"

// Get Kinetic Action
RealType Path::getK()
{
  vector<int> particles;
  for (int iP=0; iP<nPart; iP+=1)
    particles.push_back(iP);

  return getK(0,nBead,particles,0);
}

// Get Single Particle Kinetic Action
RealType Path::getK(const int iP)
{
  vector<int> particles;
  particles.push_back(iP);

  return getK(0,nBead,particles,0);
}

// Get Kinetic action for swath of path
RealType Path::getK(int b0, int b1, vector<int> &particles, int level)
{
  int nImages = 0;
  int skip = 1<<level;
  RealType levelTau = skip*tau;
  RealType tot = 0.0;
  Tvector dr(nD);
  for (int iP=0; iP<particles.size(); ++iP) {
    RealType lambda = speciesList[particles[iP]]->lambda;
    if (lambda != 0) {
      RealType i4LambdaTau = 1./(4.*lambda*levelTau);
      for (int iB=b0; iB<b1; iB+=skip) {
        Dr(bead(iP,iB),bead(iP,iB)->nextB(skip),dr);
        RealType gaussProd = 1.;
        for (int iD=0; iD<nD; iD++) {
          RealType gaussSum = 0.;
          for (int image=-nImages; image<=nImages; image++) {
            RealType dist = dr(iD) + (RealType)image*L;
            gaussSum += exp(-dist*dist*i4LambdaTau);
          }
          gaussProd *= gaussSum;
        }
        tot -= log(gaussProd);
      }
    }
  }

  return tot;
}

// Get Potential Action
RealType Path::getV()
{
  RealType Vext = 0.0;
  RealType Vint = 0.0;
  for (unsigned int iP = 0; iP < nPart-1; iP += 1) {
    for (unsigned int iB = 0; iB < nBead; iB += 1) {
      for (unsigned int jP = iP+1; jP < nPart; jP += 1) {
        Vint += getVint( bead(iP,iB) , bead(jP,iB) );
      }
      Vext += getVext( bead(iP,iB) );
    }
  }
  for (unsigned int iB = 0; iB < nBead; iB += 1) {
    Vext += getVext( bead(nPart-1,iB) );
  }

  return Vext + Vint;
}

// Get Single Particle Potential Action
RealType Path::getV( const unsigned int iP )
{
  RealType Vext = 0.0;
  RealType Vint = 0.0;
  for (unsigned int iB = 0; iB < nBead; iB += 1) {  
    for (unsigned int jP = 0; jP < iP; jP += 1)
      Vint += getVint( bead(iP,iB) , bead(jP,iB) );
    for (unsigned int jP = iP+1; jP < nPart; jP += 1) 
      Vint += getVint( bead(iP,iB) , bead(jP,iB) );

    Vext += getVext( bead(iP,iB) );
  }

  return  Vext + Vint;
}

// Get Single Bead Potential Action
RealType Path::getV( const int unsigned iP , const int iB )
{
  RealType Vext = 0.0;
  RealType Vint = 0.0;
  for (unsigned int jP = 0; jP < iP; jP += 1)
    Vint += getVint( bead(iP,iB) , bead(jP,iB) );
  for (unsigned int jP = iP+1; jP < nPart; jP += 1) 
    Vint += getVint( bead(iP,iB) , bead(jP,iB) );

  Vext += getVext( bead(iP,iB) );

  return Vext + Vint;
}

// Get Single Bead Potential Action
RealType Path::getV( Bead *bi )
{
  RealType Vext = 0.0;
  RealType Vint = 0.0;
  const unsigned int iP = bi -> p;
  for (unsigned int jP = 0; jP < iP; jP += 1)
    Vint += getVint( bi , bead(jP,bi -> b) );
  for (unsigned int jP = iP+1; jP < nPart; jP += 1) 
    Vint += getVint( bi , bead(jP,bi -> b) );

  Vext += getVext( bi );

  return Vext + Vint;
}

// Get Two Bead Interaction Action
RealType Path::getVint( Bead *b1 , Bead *b2 )
{
  return 0;
}

// Get External Potential Action (w/ LB Correction)
RealType Path::getVext( Bead *b )
{
  Tvector dr(nD);
  if(trap) {
    if(mode) {
      dr = b -> r;
    } else {
      dr = b -> rC;
    }
    return halfTauOmega2 * onePlusTau2Omega2Over12 * dot(dr, dr);
  } else {
    return 0;
  }
}

// Get Single Bead Nodal Action
RealType Path::getN( const int iP , const int iB )
{
  if (!useNodeDist) return 0;
  updateNodeDistance(bead(iP,iB));
  updateNodeDistance(bead(iP,iB)->next);
  RealType nD1 = bead(iP,iB) -> nDist;
  RealType nD2 = bead(iP,iB) -> next -> nDist;
  RealType nD1nD2 = nD1 * nD2;
  //if (!nD1) {
  //  nD1nD2 = nD2 * nD2;
  //  cerr << "nD1" << endl;
  //} else if (!nD2) {
  //  nD1nD2 = nD1 * nD1;
  //  cerr << "nD2" << endl;
  //}
  //RealType factor = -log1p(-exp(-0.5*nD1nD2*iLamTau));
  //RealType factor = -log1p(-nD1nD2*iLamTau);
  RealType factor = -log1p(-exp(-0.5*nD1nD2/(tau*bead(iP,iB)->species.lambda)));
  return factor;
}

// Get Time Slice Nodal Action
RealType Path::getNSlice( const int iB , const int skip )
{
  if (!useNodeDist) return 0;
  RealType N = 0;
  for (unsigned int iP = 0; iP < nPart; iP++) {
    N += getN(iP,iB,skip);
  }
  return N;
}

// Get Single Bead Nodal Action
RealType Path::getN( const int iP , const int iB , const int skip )
{
  if (!useNodeDist) return 0;
  Bead *b1, *b2;
  b1 = bead(iP,iB);
  b2 = bead(iP,iB)->nextB(skip);
  RealType nD1, nD2;
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
  RealType nD1nD2 = nD1 * nD2;
  if (!nD1) {
    nD1nD2 = nD2 * nD2;
  } else if (!nD2) {
    nD1nD2 = nD1 * nD1;
  }
  RealType N = 0.0;
  //N += -log1p(-exp(-nD1nD2*iLamTau/skip));
  N += -log1p(-exp(-0.5*nD1nD2/(skip*tau*bead(iP,iB)->species.lambda)));

  return N;
}

// Get Single Bead Nodal Action
RealType Path::getN( Bead *b , int skip )
{
  return getN( b->p , b->b , skip );
}

// Get Single Bead Nodal Action
RealType Path::getN( Bead *b )
{
  return getN( b -> p , b -> b );
}

// Get Nodal Action of Bead Set
RealType Path::getN( std::vector<Bead*>& beads )
{
  RealType N(0.0);
  for (std::vector<Bead*>::const_iterator b = beads.begin(); b != beads.end(); ++b)
    N += getN((*b));
  return N;
}
