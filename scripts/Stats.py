from math import sqrt
import acor

# Mean
def Mean(g):
    return sum(g)/len(g)

# Mean squared
def Mean2(g):
    return sum(g**2)/len(g)

# Sample variance
def SampVar(g):
    return Mean2(g) - Mean(g)**2

# Unbiased estimate of population variance
def Var(g):
    N = len(g)
    if N > 1:
        return ((N+0.)/(N-1.)) * SampVar(g)
    else:
        return 0.

# Standard deviation
def Sigma(g):
    return sqrt(Var(g))

# Average means
def UnweightedAvg(stats):
    mean, error, kappa = 0., 0., 0.
    for s in stats:
        mean += s[0]
        error += s[1]*s[1]
        kappa += s[2]
    N = len(stats)
    mean /= N
    error = sqrt(error)/N
    kappa /= N
    return mean, error, kappa

# Reaverage through means, errors, and N
def UnweightedReAvg(stats):
    totMean,totMean2,totError,totN = 0.,0.,0.,0.
    for s in stats:
        [mean,error,kappa,N] = s
        mean2 = error*error*(N-1.) + mean*mean
        totMean += mean*N
        totMean2 += mean2*N
        totN += N
    if totN != 0:
        totMean /= totN
        totMean2 /= totN
        var = totMean2 - totMean*totMean
        if totN > 1 and var > 0:
            totError = sqrt(var/(totN-1.))
        else:
            totError = 0.
    return totMean, totError, totN

# Statistics with auto-correlation
def stats(x):
    N = len(x)
    try:
        kappa, mean, error = acor.acor(x)
    except RuntimeError:
        kappa = 1.
        mean = Mean(x)
        var = Var(x)
        if var < 0:
            error = 0.
        else:
            error = sqrt(Var(x)/N)
    return mean, error, kappa, N
