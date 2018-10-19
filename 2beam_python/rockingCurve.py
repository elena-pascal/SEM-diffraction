import numpy as np

from TwoBeamODE import zfunc, zjac

from fileTools import readInput, writeOutput
from dynTools import thetaB, integrateOnGrid, getBackground
from strainTools import defineBetaScrew, defineBetaEdge


if __name__ == "__main__":
    # Read parameters from file.
    my_file = "rockingCurve.in"
    indata = readInput(my_file)

    thetaB_val =  thetaB(indata['g'], indata['a'], indata['c'], indata['V'])

    w = np.around(np.linspace(indata['w1'], indata['w2'], 200), decimals=2)

    # rocking curves
    bright = np.zeros(len(w))
    dark = np.zeros(len(w))
    for index, wi in enumerate(w):
        # Define inintial conditions.
        initCond = np.array([1. + 0.j, 0.+ 0.j])
        print '----------------------'
        print 'w at: ', wi

        backgr = getBackground(indata['tiltS'], \
                          indata['zmax'], indata['dstep'], \
                          zfunc, zjac, initCond, \
                          indata["Xi_0'/Xi_g'"], indata["Xi_g/Xi_g'"], wi)

        bright[index] = backgr[0]
        dark[index] = backgr[1]


    bright_RC = 'g75-3/I0_RockingCurve_' + str(indata['zmax']) + '.out'
    dark_RC = 'g75-3/Ig_RockingCurve_' + str(indata['zmax']) + '.out'


    with open(bright_RC, 'w+') as file:
        for index, wi in enumerate(w):
            file.write('%s  %s \n' % (wi, bright[index]))

    with open(dark_RC, 'w+') as file:
        for index, wi in enumerate(w):
            file.write('%s  %s \n' % (wi, dark[index]))


    print "contrast intensity written to", dark_RC
    print
