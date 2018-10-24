# Solve coupled complex ODEs using the odeintw wrapper developed by Warren Weckesser
# for the scipy odeint module which handles complex formats.
# Copyright (c) 2014, Warren Weckesser
import numpy as np
import sympy as sp

from fileTools import readInput, writeOutput, fileName
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
    print 'Bragg angle for this g is: ', thetaB_val, 'radians'
    if indata['type'] == 'Screw':
        my_beta = defineBetaScrew(indata['a'], indata['c'], indata['bs'],\
                         indata['g'], thetaB, indata['V'], indata['tiltS'], indata['rot_c'])
    elif indata['type'] == 'Edge':
        my_beta = defineBetaEdge(indata['a'], indata['c'], indata['be'],\
                         indata['rot_b'], indata['nu'], thetaB, indata['V'],\
                         indata['g'], indata['tiltS'], indata['rot_c'])
    else:
        my_beta = 0. # perfect crystal




    # Define inintial conditions.
    initCond = np.array([1. + 0.j, 0.+ 0.j])
    print '----------------------'
    print 'starting integration on grid'

    #start_time = time.time()
    # Integrate on a given grid size.
    intensityData = integrateOnGrid(indata['x_size'], indata['nx'], \
                              indata['y_size'], indata['ny'],\
                              indata['tiltS'], \
                              indata['zmax'], indata['dstep'], \
                              zfunc, zjac, initCond, \
                              indata["Xi_0/Xi_g'"], indata["Xi_g/Xi_g'"], \
                              indata["w"], my_beta)
    #end_time = time.time()
    #print integrateOnGrid.func_name, ' took', end_time - start_time, 'seconds'

    backgr = getBackground(indata['tiltS'], \
                      indata['zmax'], indata['dstep'], \
                      zfunc, zjac, initCond, \
                      indata["Xi_0/Xi_g'"], indata["Xi_g/Xi_g'"], \
                      indata["w"])

    bright = backgr[0]
    dark = backgr[1]


    # Definition of contrast
    contrast = (intensityData[0] - backgr[0]) + (intensityData[1] - backgr[1])

    contrast_file = fileName(contrast, "test")[2]

    writeOutput(contrast, contrast_file)

    print "contrast intensity written to", contrast_file
    print
