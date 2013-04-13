#include <cmath>
#include <iostream>

double expm1(double x)
{
  if (fabs(x) < 1e-5)
    return x + 0.5*x*x;
  else
    return exp(x) - 1.0;
}
