import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot

# Put in box
def PutInBox(bs,L,n_d):
    # First bead
    for d_i in range(n_d):
        n = -np.rint(bs[0][d_i]/L)
        bs[0][d_i] += n*L

    # Other beads
    for b_i in range(1,len(bs)):
        for d_i in range(n_d):
            n = -np.rint(bs[b_i][d_i]/L)
            bs[b_i][d_i] += n*L
            dist = bs[b_i][d_i] - bs[b_i-1][d_i]
            if dist > L/2.:
                bs[b_i][d_i] -= L
            elif dist < -L/2.:
                bs[b_i][d_i] += L

# Open h5 file
f = h5.File(sys.argv[1],'r')

# Get system info
nBead = np.array(f['/System/n_bead'])
L = np.array(f['/System/L'])
n_d = np.array(f['/System/n_d'])
particles = f['/System/Particles']
n_part = {}
for species in particles.keys():
    n_part[species] = np.array(particles[species+'/n_part'])

# Set up axis
fig = matplotlib.pyplot.figure()
if (n_d == 3):
    from mpl_toolkits.mplot3d import Axes3D
    ax  = fig.add_subplot(111, projection = '3d')
elif (n_d == 2):
    ax  = fig.add_subplot(111)
else:
    print 'ERROR: Can only plot 2D or 3D paths!'
    sys.exit()

# Iterate through species
colors = "bgrcmykw"
color_index = 0
path_dump = f['/Observables/PathDump']
for species in path_dump.keys():
    if species != 'type':
        n_dump = np.array(path_dump[species+'/n_dump'])
        permutation = np.array(path_dump[species+'/permutation'], dtype='int')[-1]
        positions = [p[0] for p in np.split(np.array(path_dump[species+'/positions'])[-1],n_part[species])]
        visited = []
        for i in range(n_part[species]):
            visited.append(0)
        for i in range(len(positions)):
            if not visited[i]:
                p = positions[i]
                perm_p = permutation[i][1]
                if (perm_p == i):
                    p = np.append(p, [positions[perm_p][0]], axis=0)
                    visited[i] = 1
                else:
                    while (perm_p != i):
                        p = np.append(p, positions[perm_p], axis=0)
                        perm_p = permutation[perm_p][1]
                PutInBox(p,L,n_d)
                p = np.transpose(p)
                if (n_d == 3):
                    ax.plot(p[0],p[1],zs=p[2], c=colors[color_index], marker='o', linestyle='dashed')
                elif (n_d == 2):
                    ax.plot(p[0],p[1], c=colors[color_index], marker='o', linestyle='dashed')
                else:
                    print 'ERROR: Can only plot 2D or 3D paths!'
                    sys.exit()
                color_index = (color_index+1) % len(colors)

# Show plot
ax.xaxis.set_tick_params(label1On=False)
ax.yaxis.set_tick_params(label1On=False)
if (n_d == 3):
    ax.zaxis.set_tick_params(label1On=False)
matplotlib.pyplot.show()
