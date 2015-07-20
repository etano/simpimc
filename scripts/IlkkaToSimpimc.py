import sys, os
import getopt, string
from optparse import OptionParser
import h5py as h5
import numpy as np

def cb(option, opt_str, value, parser):
    args=[]
    for arg in parser.rargs:
        if arg[0] != "-":
            args.append(arg)
        else:
            del parser.rargs[:len(args)]
        break
    if getattr(parser.values, option.dest):
        args.extend(getattr(parser.values, option.dest))
    setattr(parser.values, option.dest, args)

def main():

    # Options
    usage = "usage: python %prog [options] ../path/to/fileprefix"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--folder", action="store", type="string", dest="folder", help="Folder with Ilkka configuration files")
    parser.add_option("-o", "--out", action="store", type="string", dest="out", help="HDF5 file to store configuration")
    parser.add_option("-L", "--box", action="store", type="float", dest="L", help="Box size")
    parser.add_option("-D", "--dim", action="store", type="int", dest="D", help="Dimension")
    parser.add_option("-M", "--slices", action="store", type="int", dest="M", help="Number of time slices")
    parser.add_option("-e", "--electrons", action="callback", callback=cb, dest="electrons", default=[], help="name:num [default: %default]")
    parser.add_option("-n", "--nuclei", action="callback", callback=cb, dest="nuclei", default=[], help="name:num [default: %default]")
    (opts, args) = parser.parse_args()
    print opts,args

    NE = 0
    for e in opts.electrons:
        e = e.split(':')
        NE += int(e[1])

    NN = 0
    for n in opts.nuclei:
        n = n.split(':')
        NN += int(n[1])

    out = h5.File(opts.out,'w')
    pathDump = out.create_group('Observables').create_group('PathDump')
    pathDump.create_dataset('type',data=['PathDump'])
    i = 0
    if (NE > 0):
        positions = np.loadtxt(opts.folder+'/M'+str(opts.M)+'N'+str(NE)+'.electrons').reshape(NE,opts.M,opts.D)
        for e in opts.electrons:
            e = e.split(':')
            N = int(e[1])
            eOut = pathDump.create_group(e[0])
            eOut.create_dataset('n_dump',data=1)
            ePos0 = np.copy(positions[i:i+N])
            ePos1 = np.copy(positions[i:i+N])
            permutations = []
            for p in range(i,i+N):
                permutations.append(np.loadtxt(opts.folder+'/Perm.'+str(p+1)+'.dat', dtype=int))
            allperms = []
            permutations = np.transpose(np.array(permutations))
            for b in range(opts.M):
                for p in range(N):
                    if not (p+1 in permutations[b]):
                        print b, p+1, permutations[b]


            for p in range(N):
                print permutations[p]
                perm = permutations[0,p]-1-i
                perms = [perm]
                print p, perm
                for b in range(1,opts.M):
                    #ePos1[p,b+1] = ePos0[perm,b+1]
                    perm = permutations[b,perm]-1-i
                    perms.append(perm)
                print p, perm
                allperms.append(perms)
            sys.exit()
            eOut.create_dataset('positions',data=[ePos1])
            i += N
            
    out.flush()
    out.close()


if __name__ == "__main__":
   sys.exit(main())
