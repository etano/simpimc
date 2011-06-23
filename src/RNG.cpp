#include "RNG.h"

int seed = (int)time(0);            // random seed
CRandomSFMT1 RanGen(seed);       // make instance of random number generator
StochasticLib1 sto(seed);           // make instance of random library

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
  return RanGen.Random();
}
//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
  return (b-a)*RanGen.Random() + a;
}
//
// Generate a random integer between 1 and a given value.
// param n the largest value 
// return a uniform random value in [1,...,n]
long unifRand(long n)
{
  return RanGen.IRandom(1,n);
}
//
// Generate a box muller normal distribution random number
double normRand(double m, double s)
{				        /* mean m, standard deviation s */
	return sto.Normal(m,s);
}
////
//// Generate a box muller normal distribution random number
//double normRand(double m, double s)
//{				        /* mean m, standard deviation s */
//	float x1, x2, w, y1;
//	static float y2;
//	static int use_last = 0;
//
//	if (use_last)		        /* use value from previous call */
//	{
//		y1 = y2;
//		use_last = 0;
//	}
//	else
//	{
//		do {
//			x1 = 2.0 * RanGen.Random() - 1.0;
//			x2 = 2.0 * RanGen.Random() - 1.0;
//			w = x1 * x1 + x2 * x2;
//		} while ( w >= 1.0 );
//
//		w = sqrt( (-2.0 * log( w ) ) / w );
//		y1 = x1 * w;
//		y2 = x2 * w;
//		use_last = 1;
//	}
//
//	return( m + y1 * s );
//}
