import sys
from scipy.optimize import curve_fit
from numpy import log10, array

def usage():
    print "Usage:  %s fileName xMin xMax" % os.path.basename(sys.argv[0])

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if "-h" in argv or "--help" in argv:
        usage()
        sys.exit(2)

    # Inputs
    fileName = argv[1]
    cofactor = float(argv[2])
    xMin = float(argv[3])
    xMax = float(argv[4])
    
    # Read file
    f = open(fileName,'r')
    xs,ys = [],[]
    for line in f:
        [x,y] = [float(i) for i in line.split()]
        xs.append(x)
        ys.append(y)
    f.close()
    
    # Select data and take difference with bare Coulomb
    logxs,logys = [],[]
    for (x,y) in zip(xs,ys):
        if x > xMin and x < xMax:
            logxs.append(log10(x))
            logys.append(log10(abs(y-cofactor/x)))
    
    # Fit data
    def fn(x,m,b):
        return m*x + b
    popt,pcov = curve_fit(fn,logxs,logys)
    
    # Write fit to file
    f = open(fileName,'w')
    for (x,y) in zip(xs,ys):
        if x > xMin:
            f.write('%.10e %.10e\n' % (x, 10**fn(log10(x),*popt) + cofactor/x))
        else:
            f.write('%.10e %.10e\n' % (x,y))
    f.close()


if __name__ == "__main__":
    sys.exit(main())
