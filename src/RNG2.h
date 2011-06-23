#ifndef MYRNG_H
#define MYRNG_H
#define MULTIFILE_PROJECT

#include "StandardLibraries.h"       // Standard libraries
#include "rng/randomc.h"                   // define classes for random number generators
#include "rng/sfmt.h"                   // define classes for random number generators
#include "rng/stocc.h"                   // define classes for random number generators

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   // Include code for the chosen random number generator:
   #include "rng/sfmt.cpp"
   // define system specific user interface:
   #include "rng/userintf.cpp"
#endif

//extern int seed;            // random seed
//extern CRandomSFMT1 RanGen(seed);       // make instance of random number generator

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
extern double unifRand();
//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
extern double unifRand(double a, double b);
//
// Generate a random integer between 1 and a given value.
// param n the largest value 
// return a uniform random value in [1,...,n]
extern long unifRand(long n);
//
// Generate a box muller normal distribution random number
extern double normRand(double m, double s);

#endif
