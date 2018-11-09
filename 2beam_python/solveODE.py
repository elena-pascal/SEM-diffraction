# Solve coupled complex ODEs using the odeintw wrapper developed by Warren Weckesser
# for the scipy odeint module which handles complex formats.
# Copyright (c) 2014, Warren Weckesser
import numpy as np
import sympy as sp

from sympy.physics.vector import ReferenceFrame, express
from fileTools import readInput, writeOutput, fileName
from dynTools import thetaB, integrateOnGrid, getBackground, total_psi
from strainTools import defineBetaScrew, defineBetaEdge, defineBetaZero
from twoBeamODE import zfunc, zjac
from coordinates import sampleFrame

#--------------------------- Main starts here --------------------------------
if __name__ == "__main__":
    # Read parameters from file.
    my_file = "solveODE.in"
    indata = readInput(my_file)

    # Define static beta function.
    if indata['type'] == 'Screw':
        my_beta_num = defineBetaScrew(indata['a'], indata['c'], indata['bs'],\
                         indata['g'], thetaB, indata['V'], indata['tiltS'], indata['rot_c'])


    elif indata['type'] == 'Edge':
        my_beta_num = defineBetaEdge(indata['a'], indata['c'], indata['be'],\
                         indata['rot_b'], indata['nu'], thetaB, indata['V'],\
                         indata['g'], indata['tiltS'], indata['rot_c'])

    else:
        my_beta_num = defineBetaZero() # perfect crystal




    # Define inintial conditions.
    initCond = np.array([1. + 0.j, 0.+ 0.j])
    print '----------------------'
    print 'starting integration on grid'
    print


    # Integrate on a given grid size.
    # Return beamsArray[T[xgrid, ygrid, z], S[xgrid, ygrid, z]]
    # integrateOnGrid(x, nx, y, ny, maxZ, dZ, tiltS, func, jacob, initCond, frac0, fracg, w, beta)
    beamsArray = integrateOnGrid(indata['x_size'], indata['nx'], \
                              indata['y_size'], indata['ny'],\
                              indata['zmax'], indata['dstep'], \
                              indata['tiltS'], \
                              zfunc, zjac, initCond, \
                              indata["Xi_0/Xi_g'"], indata["Xi_g/Xi_g'"], \
                              indata["w"], my_beta_num)

    # express incident beam direction in sample frame
    LF = ReferenceFrame('LF')
    rinc_LF = -1. * LF.z
    SF = sampleFrame(indata['tiltS'], LF)
    rinc_SF = express(rinc_LF, SF)
    rinc_SF_M = np.array((rinc_SF.dot(SF.x), rinc_SF.dot(SF.y), rinc_SF.dot(SF.z)))

    g = indata['g']
    g_M = np.array((g[0], g[1], g[2]))
    # total wavefunciton on the 3D grid


    totalPsi = total_psi(beamsArray[0], beamsArray[1],   \
                        indata['x_size'], indata['nx'],  \
                        indata['y_size'], indata['ny'],  \
                        indata['zmax'], indata['dstep'], \
                        rinc_SF_M, g_M)

    # write to file totalPsi[x=0, :, :], ie x slice

    print totalPsi[0, :, :]

    realPsi = totalPsi.real
    intensityPsi = abs(totalPsi**2)
    writeOutput(realPsi[0, :, :], 'realPsiScrew.out')
    writeOutput(intensityPsi[0, :, :], 'intensityPsiScrew.out')

    print "real psi values at slice x=0 written to", 'realPsi.out'



    # background intensity beams from a perfect crystal
    backgr = getBackground(indata['tiltS'], \
                      indata['zmax'], indata['dstep'], \
                      zfunc, zjac, initCond, \
                      indata["Xi_0/Xi_g'"], indata["Xi_g/Xi_g'"], \
                      indata["w"])

    perfect_bright = backgr[0]
    perfect_dark = backgr[1]

    # total probability distribuiton defined as the quare of the wavefunction
    # is the sum of both incident squared wavefunction and diffracted squared wavefunction
    #totalI = intensityData[0] + intensityData[1]

    #integratedI_file = fileName(integratedI, "full I")[2]

    #writeOutput(totalI, integratedI_file)

    #print "integrated intensity written to", integratedI_file
    #print
