import sys
import os
import h5py as h5
import numpy as np
import matplotlib.pyplot

# Get start cut
try:
    startCut = int(sys.argv[1])
    section = sys.argv[2]
    firstArg = 3
except:
    startCut = 0
    section = sys.argv[1]
    firstArg = 2

print section

# Plot traces
fig = matplotlib.pyplot.figure()
ax  = fig.add_subplot(111)
for file in sys.argv[firstArg:]:
    f = h5.File(file,'r')
    data = np.array(f[section][startCut:])
    ax.plot(data)

# Show plot
matplotlib.pyplot.show()
