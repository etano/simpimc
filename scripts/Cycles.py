import sys
import h5py as h5
import Stats

try:
    StartCut = int(sys.argv[1])
    firstArg = 2
except:
    StartCut = 0
    firstArg = 1

stats = {}
nPart = 0
for fname in sys.argv[firstArg:]:
    f = h5.File(fname,'r')
    ys = f['Observables/Permutation/cycles'][StartCut:]
    species = f['Observables/Permutation/species'][0]
    nPart = f['System/Particles/'+species+'/nPart'][0]
    for i in range(nPart):
        if i in stats:
            stats[i].append(Stats.stats(ys[:,i]))
        else:
            stats[i] = [Stats.stats(ys[:,i])]

    f.flush()
    f.close()

for i in range(nPart):
    print i, Stats.UnweightedAvg(stats[i])
