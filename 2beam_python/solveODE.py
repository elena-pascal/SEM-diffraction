# Solve coupled complex ODEs using the odeintw wrapper developed by Warren Weckesser
# for the scipy odeint module which handles complex formats.
# Copyright (c) 2014, Warren Weckesser
import numpy as np
import sympy as sp

from fileTools import readInput, writeOutput
from dynTools import thetaB, integrateOnGrid, getBackground
from strainTools import defineBetaScrew, defineBetaEdge
from twoBeamODE import zfunc, zjac

#--------------------------- Main starts here --------------------------------
if __name__ == "__main__":
    # Read parameters from file.
    my_file = "solveODE.in"
    indata = readInput(my_file)

    # Define static beta function.
    thetaB_val =  thetaB(indata['g'], indata['a'], indata['c'], indata['V'])
    print 'thetaB is: ', thetaB_val
    if indata['type'] == 'Screw':
        my_beta = defineBetaScrew(indata['a'], indata['c'], indata['bs'],\
                         indata['g'], thetaB, indata['V'], indata['tiltS'], indata['rot_c'])
    elif indata['type'] == 'Edge':
        my_beta = defineBetaEdge(indata['a'], indata['c'], indata['be'],\
                         indata['rot_b'], indata['nu'], thetaB, indata['V'],\
                         indata['g'], indata['tiltS'], indata['rot_c'])
    else:
        my_beta = 0. # perfect crystal

    w = np.around(np.linspace(indata['w1'], indata['w2'], 100), decimals=2)

    bright = np.zeros(len(w))
    dark = np.zeros(len(w))
    for index, wi in enumerate(w):
        # Define inintial conditions.
        initCond = np.array([1. + 0.j, 0.+ 0.j])
        print '----------------------'
        print 'w at: ', wi

        #start_time = time.time()
        # Integrate on a given grid size.
        intensityData = integrateOnGrid(indata['x_size'], indata['nx'], \
                                  indata['y_size'], indata['ny'],\
                                  indata['tiltS'], \
                                  indata['zmax'], indata['dstep'], \
                                  zfunc, zjac, initCond, \
                                  indata["Xi_0/Xi_g'"], indata["Xi_g/Xi_g'"], # this fails some times and I don't know why
                                  wi, my_beta)
        #end_time = time.time()

        #print integrateOnGrid.func_name, ' took', end_time - start_time, 'seconds'


        backgr = getBackground(indata['tiltS'], \
                          indata['zmax'], indata['dstep'], \
                          zfunc, zjac, initCond, \
                          indata["Xi_0'/Xi_g'"], indata["Xi_g/Xi_g'"], wi)

        bright[index] = backgr[0]
        dark[index] = backgr[1]
        print 'depth is', backgr[2]



        # Definition of contrast
        #contrast = (intensityData[0] - getBackgound[0]) + (intensityData[1] - getBackgroung[1])


        #writeOutput(contrast, contrast_file, my_plotData[2], depth_file)
    bright_RC = 'g75-3/I0_RockingCurve_05.out'
    dark_RC = 'g75-3/Ig_RockingCurve_05.out'


    with open(bright_RC, 'w+') as file:
        for index, wi in enumerate(w):
            file.write('%s  %s \n' % (wi, bright[index]))

    with open(dark_RC, 'w+') as file:
        for index, wi in enumerate(w):
            file.write('%s  %s \n' % (wi, dark[index]))


    print "contrast intensity written to", dark_RC
    print
