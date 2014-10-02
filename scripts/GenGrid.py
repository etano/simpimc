from numpy import linspace, logspace, log10, exp

def GenGrid(o):
    rs = []
    if o['gridType']=="LINEAR":
        rs = linspace(o['r0'],o['rCut'],num=o['nGrid'],endpoint=True)
    elif o['gridType']=="LOG":
        rs = logspace(log10(o['r0']),log10(o['rCut']),num=o['nGrid'],endpoint=True)
    elif o['gridType']=="OPTIMIZED":
        rs = [o['r0']]
        f0 = 0.
        a = 10.
        dr = o['rCut']/(o['nGrid']-a)
        for iGrid in range(o['nGrid']):
            fi = f0 + 2.*(iGrid+1)*a/o['nGrid']
            rs.append(rs[iGrid] + (1. - (1./(exp(fi)+1.)))*dr)
    else:
        print 'Unrecognized grid'
    return rs
