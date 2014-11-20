import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
fig = matplotlib.pyplot.figure()
ax  = fig.add_subplot(111, projection = '3d')

# Put in box
def PutInBox(ps,L,nD):
    for p in ps:
        for i in range(nD):
            n = -np.rint(p[i]/L)
            p[i] += n*L

# Open h5 file
f = h5.File(sys.argv[1],'r')

# Get system info
nBead = np.array(f['/System/nBead'])
L = np.array(f['/System/L'])
nD = np.array(f['/System/nD'])
particles = f['/System/Particles']
nPart = {}
for species in particles.keys():
    nPart[species] = np.array(particles[species+'/nPart'])

# Iterate through species
colors = "bgrcmykw"
color_index = 0
pathDump = f['/Observables/PathDump']
for species in pathDump.keys():
    if species != 'type':
        nDump = np.array(pathDump[species+'/nDump'])
        permutation = np.array(pathDump[species+'/permutation'], dtype='int')[-1]
        positions = [p[0] for p in np.split(np.array(pathDump[species+'/positions'])[-1],nPart[species])]
        for i in range(len(positions)):
            p = positions[i]
            p = np.append(p, [positions[permutation[i][1]][0]], axis=0)
            PutInBox(p,L,nD)
            p = np.transpose(p)
            ax.plot(p[0],p[1],zs=p[2], c=colors[color_index], marker='o', linestyle='dashed')
        color_index += 1

# Show plot
matplotlib.pyplot.show()
