import sys, os
from math import sqrt

# Exact location of PAGEN scripts
PAGEN_HOME = '../../scripts/pagen'

# Units
units = {'energy':'H', 'distance':'A'}

# Constants
tau = 0.125 # Time step
L = 100.0 # Box size

# Particles
particles = [{'type': 'e', 'lambda': 0.5, 'Z': 1.0},
             {'type': 'p', 'lambda': 0.0002723089072243553, 'Z': -1.0}]

# Potential
potential = {}
potential['function'] = lambda Z1,Z2,r: Z1*Z2/r
potential['rMin'] = 0.0001 # first grid point
potential['rMax'] = 100. # last grid point
potential['nGrid'] = 1000 # number grid points
potential['gridType'] = "OPTIMIZED" # grid type (LINEAR, LOG, OPTIMIZED (Ilkka only!))

# Squarer
squarer = {}
squarer['type'] = "Ilkka" # Ilkka or David
squarer['tau'] = tau # desired timestep of PIMC simulation
squarer['D'] = 3 # dimension
squarer['rMax'] = 100.0 # maximum distance on grid
squarer['nGrid'] = 100 # number of grid points
squarer['gridType'] = "OPTIMIZED" # grid type (LINEAR, LOG, OPTIMIZED (Ilkka only!))
squarer['nSquare'] = 33 # total number of squarings to reach lowest temperature
squarer['nOrder'] = -1 # order of off-diagonal PA fit: -1 = no fit (direct spline, Ilkka only!), 0 = only diagonal, 1-3 = fit off-diagonal to 1-3 order
squarer['nTemp'] = 1 # number of temperatures for which to calculate the pair action (David only!)

# Long-range breakup
breakup = {}
breakup['D'] = 3 # dimension
breakup['L'] = L # length of box
breakup['rMin'] = 0.0001 # first grid point
breakup['rMax'] = sqrt(breakup['D'])*breakup['L']/2. # last grid point
breakup['rCut'] = breakup['L']/2. # r cutoff for ewald
breakup['nGrid'] = 1000 # number of grid points
breakup['gridType'] = "OPTIMIZED" # grid type (LINEAR, LOG, OPTIMIZED (Ilkka only!))
breakup['nKnots'] = 10 # number of knots in spline (probably fine)
breakup['nImages'] = 10 # Naive check
breakup['GenPIMCPairActionFile'] = 0 # Whether or not to generate the *.PairAction file

# Pair action objects (object type, kspace cutoff, breakup type)
# type : 2 - dU/dBeta, 1 - U, 0 - V
# breakup : 2 - Short-ranged only, 1 - Optimized breakup, 0 - Classical Ewald breakup
objects = [{'type': 0, 'breakup': 2, 'kCut': 0.},
           {'type': 1, 'breakup': 2, 'kCut': 0.},
           {'type': 2, 'breakup': 2, 'kCut': 0.}]

#-----------------------------
# DO NOT EDIT BELOW THIS LINE
#-----------------------------
sys.path.append(PAGEN_HOME)
from GenPairAction import *
run(units,particles,potential,squarer,breakup,objects)
