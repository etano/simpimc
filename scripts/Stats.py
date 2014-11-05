import sys
import os
import ctypes as ct
import numpy as np

def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def module_path():
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))

StatsLib = ct.cdll.LoadLibrary(module_path()+'/stats.so')
npct_double = np.ctypeslib.ndpointer(dtype=ct.c_double)

# Mean
StatsLib.Mean.argtypes = [npct_double, ct.c_int]
StatsLib.Mean.restype = ct.c_double
def Mean(x):
    return StatsLib.Mean(x,len(x))

# Mean squared
StatsLib.Mean2.argtypes = [npct_double, ct.c_int]
StatsLib.Mean2.restype = ct.c_double
def Mean2(x):
    return StatsLib.Mean2(x,len(x))

# Sample variance
StatsLib.Var.argtypes = [npct_double, ct.c_int]
StatsLib.Var.restype = ct.c_double
def Var(x):
    return StatsLib.Var(x,len(x))

# Standard deviation
StatsLib.StdDev.argtypes = [npct_double, ct.c_int]
StatsLib.StdDev.restype = ct.c_double
def StdDev(x):
    return StatsLib.StdDev(x,len(x))

# Average means
StatsLib.UnweightedAvg.argtypes = [npct_double, npct_double, npct_double, ct.c_int]
StatsLib.UnweightedAvg.restype = np.ctypeslib.ndpointer(dtype=ct.c_double, shape=(3,))
def UnweightedAvg(stats):
    means = stats[:,0]
    errors = stats[:,1]
    kappas = stats[:,2]
    N = len(stats)
    return StatsLib.UnweightedAvg(np.ascontiguousarray(means), np.ascontiguousarray(errors), np.ascontiguousarray(kappas), N)

# Statistics with auto-correlation
StatsLib.Stats.argtypes = [npct_double, ct.c_int]
StatsLib.Stats.restype = np.ctypeslib.ndpointer(dtype=ct.c_double, shape=(3,))
def stats(x):
    return StatsLib.Stats(x,len(x))
