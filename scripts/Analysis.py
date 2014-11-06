import sys
import os
import h5py as h5
import numpy as np
import Stats

# Return subsections
def GetSubsections(f, section_name):
    try:
        section = f[section_name]
        subsections = section.keys()
        if 'type' in subsections:
            subsections.remove('type')
        if 'data_type' in subsections:
            subsections.remove('data_type')
        return subsections
    except:
        return []

# Observable base class
class Observable:
    def __init__(self, basename, prefix, section, startCut):
        self.basename = basename
        self.prefix = prefix
        self.section = section
        self.startCut = startCut

    def Process(self, files):
        print 'Processing', self.prefix+self.section
        data = self.GetData(files)
        stats = self.DoStatistics(data)
        self.WriteToFile(stats)

    def GetData(self, files):
        print 'GetData not defined for', self.section
        return 0

    def DoStatistics(self, data):
        print 'DoStatistics not defined for', self.section
        return 0

    def WriteToFile(self, stats):
        print 'WriteToFile not defined for', self.section
        return 0

# Scalar variables
class Scalar(Observable):
    def GetData(self, files):
        data = []
        for file in files:
            f = h5.File(file,'r')
            data.append(f[self.prefix+self.section+"/x"][startCut:])
            f.flush()
            f.close()
        return data

    def DoStatistics(self, data):
        dStats = []
        for d in data:
            dStats.append(Stats.stats(d))
        stats = Stats.UnweightedAvg(np.array(dStats))
        return stats

    def WriteToFile(self, stats):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.section+'.dat','w')
        for s in stats:
            g.write('%.10e '%(s))
        g.write('\n')
        g.close()

# Histogram variables
class Histogram(Observable):
    def GetData(self, files):
        xs,ys = [],[]
        for file in files:
            f = h5.File(file,'r')
            xs = f[self.prefix+self.section+"/x"]
            ys.append(np.transpose(f[self.prefix+self.section+"/y"][startCut:]))
            f.flush()
            f.close()
        return (xs, ys)

    def DoStatistics(self, data):
        (xs, ys) = data
        stats = []
        for i in range(len(xs)):
            yStats = []
            for y in ys:
                yStats.append(Stats.stats(y[i]))
            stats.append(Stats.UnweightedAvg(np.array(yStats)))
        return (xs, stats)

    def WriteToFile(self, (xs, stats)):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.section+'.dat','w')
        for i in range(len(xs)):
            g.write(str(xs[i])+' ')
            for s in stats[i]:
                g.write('%.10e '%(s))
            g.write('\n')
        g.close()

# Pair variables
class Pair(Observable):
    def GetData(self, files):
        xs,ys = [],{}
        for file in files:
            f = h5.File(file,'r')
            xs = f[self.prefix+self.section+"/x"]
            pairs = f[self.prefix+self.section+"/y"][startCut:]
            f.flush()
            f.close()
            for pair in pairs:
                try:
                    ys[pair[0]] += pair[1]
                except:
                    ys[pair[0]] = pair[1]
        return (xs, ys)

    def DoStatistics(self, data):
        (xs, ys) = data
        tot = 0
        for x in ys.keys():
            tot += ys[x]
        for x in ys.keys():
            ys[x] = (ys[x]*1.)/(tot*1.)
        return (xs, ys)

    def WriteToFile(self, (xs, ys)):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.section+'.dat','w')
        for x in xs:
            g.write(str(x)+' ')
            try:
                g.write('%.10e \n'%(ys[x]))
            except:
                g.write('0. \n')
        g.close()

# Recursively add observables
def AddObservable(f, section, obs, startCut):

    # Check any subsections
    subsections = GetSubsections(f,section)
    for subsection in subsections:
        AddObservable(f, section+'/'+subsection, obs, startCut)

    # Try to find data type
    try:
        data_type = f[section+'/data_type'][0]
    except KeyError:
        data_type = "none"

    section_name = section.split('/')[-1]
    prefix = '/'.join(section.split('/')[:-1])
    if data_type == "scalar":
        obs.append(Scalar('./'+prefix, prefix+'/', section_name, startCut))
    elif data_type == "histogram":
        obs.append(Histogram('./'+prefix, prefix+'/', section_name, startCut))
    elif data_type == "pairs":
        obs.append(Pair('./'+prefix, prefix+'/', section_name, startCut))

# Get start cut
try:
    startCut = int(sys.argv[1])
    firstArg = 2
except:
    startCut = 0
    firstArg = 1

# Get h5 files
files = []
for file in sys.argv[firstArg:]:
    files.append(file)

# Get observables
obs = []
f0 = h5.File(files[0],'r')
ob_names = GetSubsections(f0,'/')
for ob_name in ob_names:
    AddObservable(f0, ob_name, obs, startCut)
f0.flush()
f0.close()

# Compute statistics
for ob in obs:
    ob.Process(files)
