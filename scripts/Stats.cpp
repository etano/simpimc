#include <math.h>
#include <iostream>

extern "C" {

  double Mean(double* x, int N)
  {
    double tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i];
    return tot/N;
  }

  double Mean2(double* x, int N)
  {
    double tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i]*x[i];
    return tot/N;
  }

  double RootMean2(double* x, int N)
  {
    double tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i]*x[i];
    return sqrt(tot)/N;
  }

  double Var(double* x, int N)
  {
    double m2 = Mean2(x,N);
    double m = Mean(x,N);
    return m2 - m*m;
  }

  double StdDev(double* x, int N)
  {
    return sqrt(Var(x,N));
  }

  double C(double* x, int t, double mean, double var, int N)
  {
    if (!var)
      return 0.;
    double tot = 0.;
    for (int i=0; i<N-t; ++i)
      tot += (x[i]-mean)*(x[N-t-i]-mean);
    return tot/(var*(N-t));
  }

  double Kappa(double* x, double mean, double var, int N)
  {
    double tot = 0.;
    for (int i=1; i<N; ++i) {
      double c = C(x,i,mean,var,N);
      if (c <= 0)
        break;
      else
        tot += c;
    }
    return 1. + 2.*tot;
  }

  double* UnweightedAvg(double* means, double* errors, double* kappas, int N)
  {
    double* stats = new double[3];
    stats[0] = Mean(means,N);
    stats[1] = RootMean2(errors,N);
    stats[2] = Mean(kappas,N);
    return stats;
  }

  double* Stats(double* x, int N)
  {
    double mean = Mean(x,N);
    double var = Var(x,N);
    double kappa = Kappa(x, mean, var, N);
    int Neff = N/kappa;
    double error = sqrt(var/Neff);

    double* stats = new double[3];
    stats[0] = mean;
    stats[1] = error;
    stats[2] = kappa;
    return stats;
  }

}
