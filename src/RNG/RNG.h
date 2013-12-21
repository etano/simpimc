#ifndef RNGClass_H
#define RNGClass_H

#include "../config.h"

#define MULTIFILE_PROJECT

#include "randomc.h"                   // define classes for random number generators
#include "sfmt.h"                   // define classes for random number generators
#include "stocc.h"                   // define classes for random number generators

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   // Include code for the chosen random number generator:
   #include "sfmt.cpp"
   // define system specific user interface:
   #include "userintf.cpp"
#endif

class RNG
{
private:
  CRandomSFMT1 RanGen;  // make instance of random number generator
  StochasticLib1 sto;  // make instance of random library
protected:
public:
  // Constructor
  RNG(int seed)
    : RanGen(seed), sto(seed)
  {
    std::cout << "RNG seed: " << seed << "\n";
  }

  // Random Functions
  double unifRand();
  double unifRand( const double a , const double b );
  long unifRand( const long n );
  void unifRand( Tvector& r );
  void unifRand( Tvector& r , const double l );
  double normRand( const double m , const double s );
  void normRand( Tvector& r , const double m , const double s );
};

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
inline double RNG::unifRand()
{
  return RanGen.Random();
}

// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
inline double RNG::unifRand( const double a , const double b )
{
  return (b-a)*RanGen.Random() + a;
}

// Generate a random integer between 1 and a given value.
// param n the largest value 
// return a uniform random value in [1,...,n]
inline long RNG::unifRand( const long n )
{
  return RanGen.IRandom(1,n);
}

// Generate a uniform random vector of length 1
inline void RNG::unifRand( Tvector& r )
{
  for (unsigned int iD = 0; iD < r.n_elem; iD += 1) 
    r(iD) = unifRand(-1,1);
  r /= norm(r, 2);
}

// Generate a uniform random vector of length l
inline void RNG::unifRand( Tvector& r , const double l )
{
  unifRand(r);
  r *= l;
}

// Generate a normal distribution random number
inline double RNG::normRand( const double m , const double s )
{				        /* mean m, standard deviation s */
	return sto.Normal(m,s);
}

// Generate a normal random vector of length 1
inline void RNG::normRand( Tvector& r , const double m , const double s )
{
  for (unsigned int iD = 0; iD < r.n_elem; iD += 1) 
    r(iD) = normRand(m,s);
}

#endif
