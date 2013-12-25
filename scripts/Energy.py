import sys
import numpy as np
import h5py as h5
import Stats

StartCut = 1

ENames = ['Kinetic']
EStats = {}
for EName in ENames:
  EStats[EName] = []
for fname in sys.argv[1:]:
  print fname
  f = h5.File(fname,'r')

  for EName in ENames:
    try:
      Es = f['Observables/Energy/'+EName][StartCut:]
      EStats[EName].append(Stats.stats(Es))
    except:
      pass

  f.flush()
  f.close()


for EName in ENames:
  try:
    TotStats = Stats.UnweightedAvg(EStats[EName])
    print EName, TotStats
  except:
    pass

