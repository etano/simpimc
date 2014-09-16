import sys
import h5py as h5
import Stats

try:
    StartCut = int(sys.argv[1])
    firstArg = 2
except:
    StartCut = 0
    firstArg = 1

EStats = {}
for fname in sys.argv[firstArg:]:
    f = h5.File(fname,'r')
    ENames = f['Observables/Energy'].keys()
    for EName in ENames:
        try:
            Es = f['Observables/Energy/'+EName][StartCut:]
            if EName in EStats:
                EStats[EName].append(Stats.stats(Es))
            else:
                EStats[EName] = [Stats.stats(Es)]
        except:
            pass

    f.flush()
    f.close()

for EName in EStats.keys():
    try:
        TotStats = Stats.UnweightedAvg(EStats[EName])
        print EName, TotStats
    except:
        pass

