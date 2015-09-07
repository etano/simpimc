import sys
import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

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
fig = plt.figure()
ax  = fig.add_subplot(111)
for file in sys.argv[firstArg:]:
    f = h5.File(file,'r')
    data = np.array(f[section][startCut:])
    ax.plot(data,label=file)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show plot
plt.title(section)
plt.savefig('fig.eps')
