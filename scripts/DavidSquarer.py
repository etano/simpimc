import sys
import os
import subprocess
from math import sqrt, pi
from numpy import loadtxt
import DavidParse
from GenGrid import GenGrid

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

def GenPotgenInput(prefix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,r0,rCut,kCut,nTemp,nSquare,breakup):
    # Determine temperatures
    if nTemp > 8:
        print 'WARNING: Max nTemp is 8!'
        nTemp = 8
    minTau = tau
    maxTau = minTau*(2**(nTemp-1))
    if nSquare < 14 + nTemp-1:
        nSquare = 14 + nTemp-1

    print 'Creating '+prefix+'.in'
    f = open(prefix+'.in','w')
    f.write('UNITS '+unit1+' '+unit2+'\n')
    f.write('TYPE '+type1+' %f\n' % (lam1))
    f.write('TYPE '+type2+' %f\n' % (lam2))
    f.write('GRID %i %s %f %f\n' % (nGrid,gridType,r0,rCut))
    f.write('SQUARER %f %i %i 3 30 %i\n' % (1./maxTau,nTemp,D,nSquare))
    boxString = ' '.join([str(L) for i in range(D)])
    if breakup == 2:
        f.write('POT COUL %f %f 0.D0\n' % (rCut,Z1Z2))
    if breakup == 1:
        f.write('POT COULOPT %f %f %f %i 0.D0 1.0 %s\n' % (rCut,kCut,Z1Z2,D,boxString))
    elif breakup == 0:
        f.write('POT COULLR %f %f %f %i 0.D0 1.0 %s\n' % (rCut,kCut,Z1Z2,D,boxString))
    f.close()

def GenPairActionInput(prefix,type1,lam1,type2,lam2,D,longRange):
    f = open(prefix+'.PairAction','w')
    f.write('Section (Fits)')
    f.write('\n{')
    f.write('\n  string Type="DavidFit";')
    f.write('\n  int NumOffDiagonalTerms = 3;')
    f.write('\n  Section (Particle1)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type1+'";')
    f.write('\n    double lambda ='+str(lam1)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  Section (Particle2)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type2+'";')
    f.write('\n    double lambda ='+str(lam2)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  string Daviddmfile = "'+prefix+'.h5";')
    if longRange:
        f.write('\n  bool longRange = true;')
    f.write('\n}')
    f.close()

def Square(particles,squarer,objects):
    for i in xrange(0, len(particles)):
        for j in xrange(i, len(particles)):

            # Assign particle attributes
            [type1, lam1, Z1] = particles[i]['type'], particles[i]['lambda'], particles[i]['Z']
            [type2, lam2, Z2] = particles[j]['type'], particles[j]['lambda'], particles[j]['Z']
            print '****************************************'
            print '****', type1, ', lam1 =', lam1, ', Z1 =', Z1
            print '****', type2, ', lam2 =', lam2, ', Z2 =', Z2

            prefix = type1+'-'+type2+'.sq'

            # Squarer
            print 'Performing squaring procedure...'
            #subprocess.call(['squarer',prefix])

            # Density Matrix Parser
            print 'Parsing density matrix...'
            #DavidParse.main(['',prefix+'.dm'])

            # Build PairAction File
            print 'Creating pairAction file...'
            longRange = 0
            if objects[0]['breakup'] != 2:
                longRange = 1
            GenPairActionInput(prefix,type1,lam1,type2,lam2,squarer['D'],longRange)

def Breakup(units,particles,potential,squarer,breakup,objects):
    paIndex = 0
    for i in xrange(0, len(particles)):
        for j in xrange(i, len(particles)):
            # Start from 1
            paIndex += 1

            # Assign particle attributes
            [type1, lam1, Z1] = particles[i]['type'], particles[i]['lambda'], particles[i]['Z']
            [type2, lam2, Z2] = particles[j]['type'], particles[j]['lambda'], particles[j]['Z']
            print '****************************************'
            print '****', type1, ', lam1 =', lam1, ', Z1 =', Z1
            print '****', type2, ', lam2 =', lam2, ', Z2 =', Z2

            prefix = type1+'-'+type2+'.sq'

            # Generate potgen input
            GenPotgenInput(prefix,units['energy'],units['distance'],type1,type2,lam1,lam2,Z1*Z2,
                           breakup['L'],breakup['D'],squarer['tau'],breakup['gridType'],
                           breakup['nGrid'],breakup['r0'],breakup['rCut'],objects[0]['kCut'],
                           squarer['nTemp'],squarer['nSquare'],objects[0]['breakup'])

            # Write potential
            f = open('v.'+str(paIndex)+'.txt','w')
            rs = GenGrid(potential)
            for r in rs:
                f.write('%.10E %.10E\n' % (r,potential['function'](Z1,Z2,r)))
            f.close()

            # Write grid
            f = open('grid.'+str(paIndex)+'.txt','w')
            rs = GenGrid(breakup)
            for r in rs:
                f.write('%.10E\n' % (r))
            f.close()

            # Do breakup
            if objects[0]['breakup'] != 2:
                if breakup['gridType']=="LINEAR":
                    gridIndex = 0
                elif breakup['gridType']=="LOG":
                    gridIndex = 1
                elif breakup['gridType']=="OPTIMIZED":
                    gridIndex = 2
                else:
                    print 'Unrecognized grid:', breakup['gridType']

                subprocess.call(['ewald',str(breakup['L']),str(objects[0]['kCut']),str(breakup['r0']),
                                 str(breakup['rCut']),str(breakup['nGrid']),str(gridIndex),
                                 str(Z1*Z2),str(objects[0]['breakup']),str(objects[0]['type']),str(paIndex),
                                 str(breakup['nKnots']),str(squarer['tau']),str(breakup['nImages'])])

                # Write .yk file
                kData = loadtxt('v.'+str(paIndex)+'.k.txt')
                ks = kData[:,0]
                Vks = kData[:,1]
                f = open(prefix+'.yk','w')
                for [k,Vk] in zip(ks,Vks):
                    f.write('  %.10E'%k+'       %.10E'%Vk+'\n')
                f.close()
            else:
                f = open('v.'+str(paIndex)+'.r.txt','w')
                rs = GenGrid(potential)
                for r in rs:
                    f.write('%.10E %.10E %.10E\n' % (r,potential['function'](Z1,Z2,r),0.))
                f.close()

            # Write .dm and .in file for squarer
            g = open(prefix+'.dm','w')
            f = open(prefix+'.in','r')
            for line in f:
                g.write(line)
            f.close()
            if objects[0]['breakup'] != 2:
                g.write('\n VIMAGE ')
                if (breakup['D'] == 2):
                    g.write(str(-3.90026492*Z1*Z2/breakup['L']))
                elif (breakup['D'] == 3):
                    g.write(str(-2.837297479*Z1*Z2/breakup['L']))
            g.write('\n POTTAIL 0.0')
            g.write('\n RANK 2 '+str(breakup['nGrid'])+' 1')
            g.write('\n GRID 1 '+breakup['gridType']+' '+str(breakup['r0'])+' '+str(breakup['rCut']))
            g.write('\n LABEL 1 r')
            g.write('\n BEGIN potential 0')
            rData = loadtxt('v.'+str(paIndex)+'.r.txt')
            count = 0
            g.write('\n  ')
            for [r,V] in rData:
                if r != 0.:
                    if objects[0]['breakup'] != 2:
                        g.write('%.10E'%(potential['function'](Z1,Z2,r) - VLong)+'  ')
                    else:
                        g.write('%.10E'%(V)+'  ')
                    count += 1
                    if count % 5 == 0:
                        g.write('\n  ')
            g.write('\n')
            g.close()
