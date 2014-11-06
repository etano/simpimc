#include <math.h>

typedef double RealType;
extern "C" {

  RealType Mean(RealType* x, int N)
  {
    RealType tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i];
    return tot/N;
  }

  RealType Mean2(RealType* x, int N)
  {
    RealType tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i]*x[i];
    return tot/N;
  }

  RealType RootMean2(RealType* x, int N)
  {
    RealType tot = 0.;
    for (int i=0; i<N; ++i)
      tot += x[i]*x[i];
    return sqrt(tot)/N;
  }

  RealType Var(RealType* x, int N)
  {
    RealType m2 = Mean2(x,N);
    RealType m = Mean(x,N);
    return m2 - m*m;
  }

  RealType StdDev(RealType* x, int N)
  {
    return sqrt(Var(x,N));
  }

  RealType C(RealType* x, int t, RealType mean, RealType var, int N)
  {
    if (!var)
      return 0.;
    RealType tot = 0.;
    for (int i=0; i<N-t; ++i)
      tot += (x[i]-mean)*(x[N-t-i]-mean);
    return tot/(var*(N-t));
  }

  RealType Kappa(RealType* x, RealType mean, RealType var, int N)
  {
    RealType tot = 0.;
    for (int i=1; i<N; ++i) {
      RealType c = C(x,i,mean,var,N);
      if (c <= 0)
        break;
      else
        tot += c;
    }
    return 1. + 2.*tot;
  }

  RealType* UnweightedAvg(RealType* means, RealType* errors, RealType* kappas, int N)
  {
    RealType* stats = new RealType[3];
    stats[0] = Mean(means,N);
    stats[1] = RootMean2(errors,N);
    stats[2] = Mean(kappas,N);
    return stats;
  }

  RealType* Stats(RealType* x, int N)
  {
    RealType mean = Mean(x,N);
    RealType var = Var(x,N);
    RealType kappa = Kappa(x, mean, var, N);
    int Neff = N/kappa;
    RealType error = sqrt(var/Neff);

    RealType* stats = new RealType[3];
    stats[0] = mean;
    stats[1] = error;
    stats[2] = kappa;
    return stats;
  }

}
