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
    def __init__(self, prefix, name, type, startCut):
        self.basename = './'+prefix
        self.prefix = prefix+'/'
        self.name = name
        self.type = type
        self.startCut = startCut

    def Process(self, files):
        print 'Processing', self.prefix+self.name
        stats = self.GetDataStats(files)
        self.WriteToFile(stats)

    def GetDataStats(self, files):
        print 'GetDataStats not defined for', self.name
        return 0

    def WriteToFile(self, stats):
        print 'WriteToFile not defined for', self.name
        return 0

    def AdjustBySign(self, sign_data):
        print 'AdjustBySign not defined for', self.name
        return 0

# Scalar variables
class Scalar(Observable):
    def GetDataStats(self, files):
        data = []
        for file in files:
            f = h5.File(file,'r')
            data.append(Stats.stats(np.array(f[self.prefix+self.name+"/x"][self.startCut:])))
            f.flush()
            f.close()
        stats = Stats.UnweightedAvg(np.array(data))
        return stats

    def WriteToFile(self, stats):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.name+'.dat','w')
        for s in stats:
            g.write('%.10e '%(s))
        g.write('\n')
        g.close()

    def AdjustBySign(self, sign_data):
        data = np.loadtxt(self.basename+'/'+self.name+'.dat')
        mean = data[0]/sign_data[0]
        if data[0] != 0:
            err = mean*np.sqrt(pow(data[1]/data[0],2) + pow(sign_data[1]/sign_data[0],2))
        else:
            err = 0.
        kappa = data[2]
        g = open(self.basename+'/'+self.name+'.adj.dat','w')
        g.write('%.10e %.10e %.10e\n'%(mean,err,kappa))
        g.close()

# Histogram variables
class Histogram(Observable):
    def GetDataStats(self, files):
        xs,yStats = [],[]
        count = 0
        for j in range(len(files)):
            f = h5.File(files[j],'r')
            xs = np.array(f[self.prefix+self.name+"/x"])
            ys = np.transpose(f[self.prefix+self.name+"/y"][self.startCut:])
            f.flush()
            f.close()
            yStats.append([])
            for i in range(len(xs)):
                yStats[j].append(Stats.stats(np.array(ys[i])))
        stats = []
        for i in range(len(xs)):
            yStatsi = [x[i] for x in yStats]
            stats.append(Stats.UnweightedAvg(np.array(yStatsi)))
        return (xs, stats)

    def WriteToFile(self, (xs, stats)):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.name+'.dat','w')
        for i in range(len(xs)):
            g.write(str(xs[i])+' ')
            for s in stats[i]:
                g.write('%.10e '%(s))
            g.write('\n')
        g.close()

    def AdjustBySign(self, sign_data):
        data = np.loadtxt(self.basename+'/'+self.name+'.dat')
        xs = data[:,0]
        means = data[:,1]/sign_data[0]
        errs = data[:,2]
        for i in range(len(errs)):
            if data[i,1] != 0:
                errs[i] = np.abs(means[i])*np.sqrt(pow(errs[i]/data[i,1],2) + pow(sign_data[1]/sign_data[0],2))
        kappas = data[:,3]
        g = open(self.basename+'/'+self.name+'.adj.dat','w')
        for (x,mean,err,kappa) in zip(xs,means,errs,kappas):
            g.write('%.10e %.10e %.10e %.10e\n'%(x,mean,err,kappa))
        g.close()

# Pair variables
class Pair(Observable):
    def GetDataStats(self, files):
        xs,ys = [],{}
        for file in files:
            f = h5.File(file,'r')
            xs = np.array(f[self.prefix+self.name+"/x"])
            pairs = np.array(f[self.prefix+self.name+"/y"][self.startCut:])
            f.flush()
            f.close()
            for pair in pairs:
                try:
                    ys[pair[0]] += pair[1]
                except:
                    ys[pair[0]] = pair[1]
        tot = 0
        for x in ys.keys():
            tot += ys[x]
        for x in ys.keys():
            ys[x] = (ys[x]*1.)/(tot*1.)
        return (xs, ys)

    def WriteToFile(self, (xs, ys)):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.name+'.dat','w')
        for x in xs:
            g.write(str(x)+' ')
            try:
                g.write('%.10e \n'%(ys[x]))
            except:
                g.write('0. \n')
        g.close()

# Avergaed pair variables
class AvgPair(Observable):
    def GetDataStats(self, files):
        xs,ys = [],{}
        for file in files:
            f = h5.File(file,'r')
            xs = np.array(f[self.prefix+self.name+"/x"])
            avg_pairs = np.array(f[self.prefix+self.name+"/y"][self.startCut:])
            f.flush()
            f.close()
            for avg_pair in avg_pairs:
                sector = avg_pair[0]
                [eM,varM,M] = avg_pair[1:]
                try:
                    [eN,varN,N] = ys[sector]
                    NM = N + M
                    eNM = (N*eN + M*eM)/NM
                    varNM = ((N*varN + N*eN*eN + M*varM + M*eM*eM)/NM) - eNM*eNM
                    ys[sector] = [eNM,varNM,NM]
                except:
                    ys[sector] = [eM,varM,M]
        for x in ys:
            ys[x][1] = np.sqrt(ys[x][1]/ys[x][2])
        return (xs, ys)

    def WriteToFile(self, (xs, ys)):
        if not os.path.exists(self.basename):
            os.makedirs(self.basename)
        g = open(self.basename+'/'+self.name+'.dat','w')
        for x in xs:
            g.write(str(x)+' ')
            try:
                g.write('%.10e %.10e %i \n'%(ys[x][0], ys[x][1], ys[x][2]))
            except:
                g.write('0. 0. 0 \n')
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

    # Try to find observable type
    try:
        section_type = f[section+'/type'][0]
    except KeyError:
        section_type = "none"

    section_name = section.split('/')[-1]
    prefix = '/'.join(section.split('/')[:-1])
    if data_type == "scalar":
        obs.append(Scalar(prefix, section_name, section_type, startCut))
    elif data_type == "histogram":
        obs.append(Histogram(prefix, section_name, section_type, startCut))
    elif data_type == "pairs":
        obs.append(Pair(prefix, section_name, section_type, startCut))
    elif data_type == "avg_pairs":
        obs.append(AvgPair(prefix, section_name, section_type, startCut))

def usage():
    print "Usage: %s [start cut] file0.h5 file1.h5 ..." % os.path.basename(sys.argv[0])
    sys.exit(2)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if "-h" in argv or "--help" in argv:
        usage()

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
    adjustBySign = 0
    sign_data = np.array((1,3))
    for ob in obs:
        ob.Process(files)
        if ob.type == 'Sign':
            adjustBySign = 1
            sign_data = np.loadtxt(ob.basename+'/'+ob.name+'.dat')

    # Adjust by sign
    if adjustBySign:
        for ob in obs:
            if ('Observables' in ob.basename) and (not ('Time' in ob.basename)) and (not ('Permutation' in ob.basename)) and (ob.type != 'Sign'):
                print 'Adjusting', ob.name, 'by sign weight'
                ob.AdjustBySign(sign_data)

if __name__ == "__main__":
  sys.exit(main())
