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
for fname in sys.argv[firstArg:]:
    f = h5.File(fname,'r')
    try:
        grs = f['Observables/PC_ep/gr'][StartCut:][:,:,0]
        r = f['Observables/PC_ep/r'][:,0]
        for i in range(len(r)):
            if i in stats:
                stats[i].append(Stats.stats(grs[:,i]))
            else:
                stats[i] = [Stats.stats(grs[:,i])]
    except:
        pass

    f.flush()
    f.close()

for i in range(len(r)):
    try:
        print r[i], Stats.UnweightedAvg(stats[i])
    except:
        pass
