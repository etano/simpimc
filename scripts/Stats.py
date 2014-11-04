from math import sqrt
import acor

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
    kappa, mean, error = acor.acor(x)
    return mean, error, kappa, N
