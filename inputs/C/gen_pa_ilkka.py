import sys, os
from math import sqrt

# Exact location of PAGEN scripts
PAGEN_HOME = '../../scripts/pagen'
sys.path.append(PAGEN_HOME)
from GenPairAction import *

# Settings
D = 3 #dimension
NE = 20
NP = 20
NC = 1
T_Kelvin = 10.0e6 # temperature in K
mass_density = 100.0 # g/cm^3
M = 16 # number of time slices

# Calculate things
KelvinPerHartree = 315774.64
T = T_Kelvin/KelvinPerHartree # temperature in Hartree
tau = (1/T)/M # used tau
E_mass = 9.10938291e-28 # g
amu = 1822.888486192*E_mass # g
P_mass = 1.007825*amu # g
C_mass = 12*amu # g
a0cm = 5.2917721092e-9 # cm
L = pow((NE*E_mass + NP*P_mass + NC*C_mass)/mass_density,1./3.)/a0cm
kCut = 14./(L/2.) # k-space cutoff
print 'tau', tau, 'kCut', kCut, 'L', L, 'beta', M*tau

# Species
e = {'type': 'e', 'lambda': 0.5, 'Z': -1.0}
p = {'type': 'p', 'lambda': 0.5*E_mass/P_mass, 'Z': 1.0}
C = {'type': 'C', 'lambda': 0.5*E_mass/C_mass, 'Z': 12.0}
print e, p, C

# Potential
potential = {}
potential['function'] = lambda Z1,Z2,r: Z1*Z2/r
potential['r_min'] = 0.0001 # first grid point
potential['r_max'] = 100. # last grid point
potential['n_grid'] = 10000 # number grid points
potential['grid_type'] = "OPTIMIZED" # grid type (LINEAR, LOG, LOGLIN (David only!), OPTIMIZED (Ilkka only!))

# Squarer
squarer = {}
squarer['type'] = "Ilkka" # Ilkka, David, or None
squarer['tau'] = tau # desired timestep of PIMC simulation
squarer['n_d'] = D # dimension
squarer['r_max'] = 100.0 # maximum distance on grid
squarer['n_grid'] = 100 # number of grid points
squarer['grid_type'] = "OPTIMIZED" # grid type (LINEAR, LOG, LOGLIN (David only!), OPTIMIZED (Ilkka only!))
squarer['n_square'] = 33 # total number of squarings to reach lowest temperature
squarer['n_order'] = -1 # order of off-diagonal PA fit: -1 = no fit (direct spline, Ilkka only!), 0 = only diagonal, 1-3 = fit off-diagonal to 1-3 order
squarer['n_temp'] = 1 # number of temperatures for which to calculate the pair action (David only!)

# Squarer
bare_squarer = {}
bare_squarer['type'] = "None" # Ilkka, David, or None
bare_squarer['tau'] = tau # desired timestep of PIMC simulation
bare_squarer['n_d'] = D # dimension
bare_squarer['r_max'] = 100.0 # maximum distance on grid
bare_squarer['n_grid'] = 100 # number of grid points
bare_squarer['grid_type'] = "OPTIMIZED" # grid type (LINEAR, LOG, LOGLIN (David only!), OPTIMIZED (Ilkka only!))

# Long-range breakup
breakup = {}
breakup['type'] = 'StandardEwald' # StandardEwald, StandardEwald, or None
breakup['n_d'] = D # dimension
breakup['L'] = L # length of box
breakup['tau'] = tau # desired timestep of PIMC simulation
breakup['r_min'] = 0.0001 # first grid point
breakup['r_max'] = sqrt(breakup['n_d'])*breakup['L']/2. # last grid point
breakup['r_cut'] = breakup['L']/2. # r cutoff for ewald
breakup['k_cut'] = 14./(L/2.) # k cutoff for ewald
breakup['n_grid'] = 1000 # number of grid points
breakup['grid_type'] = "OPTIMIZED" # grid type (LINEAR, LOG, LOGLIN (David only!), OPTIMIZED (Ilkka only!))
breakup['n_knots'] = 10 # number of knots in spline (probably fine)
breakup['n_images'] = 10 # Naive check

# Pair action objects
pa_objects = [
{'species_a': e, 'species_b': e, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': e, 'species_b': p, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': e, 'species_b': C, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': p, 'species_b': p, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': p, 'species_b': C, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': C, 'species_b': C, 'potential': potential, 'breakup': breakup, 'squarer': bare_squarer}
]

# Run
run(pa_objects)
