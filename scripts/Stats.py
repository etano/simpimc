import sys
import os
import ctypes
import numpy as np

def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def module_path():
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))

StatsLib = ctypes.cdll.LoadLibrary(module_path()+'/stats.so')

# Mean
StatsLib.Mean.argtypes = [ctypes.c_void_p,ctypes.c_int]
StatsLib.Mean.restype = ctypes.c_double
def Mean(x):
    return StatsLib.Mean(x.ctypes,len(x))

# Mean squared
StatsLib.Mean2.argtypes = [ctypes.c_void_p,ctypes.c_int]
StatsLib.Mean2.restype = ctypes.c_double
def Mean2(x):
    return StatsLib.Mean2(x.ctypes,len(x))

# Sample variance
StatsLib.Var.argtypes = [ctypes.c_void_p,ctypes.c_int]
StatsLib.Var.restype = ctypes.c_double
def Var(x):
    return StatsLib.Var(x.ctypes,len(x))

# Standard deviation
StatsLib.StdDev.argtypes = [ctypes.c_void_p,ctypes.c_int]
StatsLib.StdDev.restype = ctypes.c_double
def StdDev(x):
    return StatsLib.StdDev(x.ctypes,len(x))

# Average means
StatsLib.UnweightedAvg.argtypes = [ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int]
StatsLib.UnweightedAvg.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(3,))
def UnweightedAvg(stats):
    means = stats[:,0]
    errors = stats[:,1]
    kappas = stats[:,2]
    N = len(stats)
    return StatsLib.UnweightedAvg(means.ctypes,errors.ctypes,kappas.ctypes,N)

# Statistics with auto-correlation
StatsLib.Stats.argtypes = [ctypes.c_void_p,ctypes.c_int]
StatsLib.Stats.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(3,))
def stats(x):
    return StatsLib.Stats(x.ctypes,len(x))
