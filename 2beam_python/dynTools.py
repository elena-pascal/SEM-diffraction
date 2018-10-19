import sympy as sp
import numpy as np

#from odeintw import odeintw
import odeintw_cython# import odeintw
import time


#--------------------------- Functions start here ----------------------------



#@profile
def integrateOnGrid(x, nx, y, ny, tiltS, maxZ, dZ, func, jacob, initCond, frac0, fracg, w, beta):
    '''
    Integrate a set of complex coupled ODEs on a given grid.
    Returns bright field intensity, dark field intensity and penetration depth pixel binned.
    '''
    xM = np.linspace(-x*0.5, x*0.5, int(nx))
    # tilt correction
    yM = np.linspace(-y*0.5/np.cos(np.radians(tiltS)), y*0.5/np.cos(np.radians(tiltS)), int(ny))


    SumIg = np.zeros((int(nx), int(ny)))
    SumI0 = np.zeros((int(nx), int(ny)))
    Depth = np.zeros((int(nx), int(ny)))

    for xidx, xgrid in np.ndenumerate(xM): # scan across x
        print "x at", xgrid
        for yidx, ygrid in np.ndenumerate(yM): # scan across y
            # 2 beam initial conditions.
            phi0 = initCond

            reachedMaxDepth = False
            tiltCorrect = ygrid * np.tan(np.radians(tiltS))
            depth = 0. - tiltCorrect

            while (not reachedMaxDepth) and (depth < (maxZ - tiltCorrect)): # integrate intensity along a column
                # Calculate beta component at this position and depth.
                localBeta = beta(xgrid, ygrid, depth)
                # From known initial condition at the start of this slice
                # calculate the beam intensity at the exit of the slice of width dstep
                t = np.linspace(depth , depth+dZ, 101)

                #starttime = time.time()
                # Call odeintw.
                phi, infodict = odeintw_cython.odeintw(func, phi0, t, args=(frac0, fracg, w, localBeta),
                                Dfun=jacob, full_output=True)
                #endtime = time.time()
                #print 'odeintw took', endtime - starttime, 'seconds'

                absI0 =  abs(phi[:,0])
                absIg =  abs(phi[:,1])

                # Add intensity to beam sums.
                SumI0[xidx, yidx] += (absI0*absI0).sum()
                SumIg[xidx, yidx] += (absIg*absIg).sum()

                #totalI = abs(endPhi[0])*abs(endPhi[0]) + abs(endPhi[1])*abs(endPhi[1])
                totalI = absI0[100]*absI0[100] + absIg[100]*absIg[100]


                # If more than 99% of initial intensity is lost
                # consider it reached max penetration depth
                if totalI < 0.01:
                    #print "max depth", depth + tiltCorrect
                    reachedMaxDepth = True

                    # Save final, tilt-corrected, depth values.
                    Depth[xidx, yidx] = (depth + tiltCorrect) * np.sin(-np.radians(tiltS))

                # Reinitialise the initial condition with calulated one
                phi0 = np.array([phi[100,0], phi[100,1]])

                # Move one step lower
                depth += dZ
    return [SumI0, SumIg, Depth]



def getBackground(tiltS, maxZ, dZ, func, jacob, initCond, frac0, fracg, w):
    ''' Calculate intensities for a perfect crystal for backgound substraction.'''
    print "Calculating background..."
    # 2 beam initial conditions.
    phi0 = initCond

    reachedMaxDepth = False
    ygrid = 0.
    tiltCorrect = ygrid * np.tan(np.radians(tiltS))
    depth = 0. - tiltCorrect

    SumI0_backgr = 0.
    SumIg_backgr = 0.

    Beta_backgr = 0.
    #while (not reachedMaxDepth) and (depth < (maxZ - tiltCorrect)): # integrate intensity along a column
    while  (depth < (maxZ - tiltCorrect)):
        # From known initial condition at the start of this slice
        # calculate the beam intensity at the exit of the slice of width dstep
        tfar = np.linspace(depth , depth+dZ, 101)

        phi_backgr, infodict = odeintw_cython.odeintw(func, phi0, tfar, args=(frac0, fracg, w, Beta_backgr),
                        Dfun=jacob, full_output=True)

        absI0_backgr =  abs(phi_backgr[:,0])
        absIg_backgr =  abs(phi_backgr[:,1])

        #   Add intensity to beam sums.
        SumI0_backgr += (absI0_backgr*absI0_backgr).sum()
        SumIg_backgr += (absIg_backgr*absIg_backgr).sum()

        #totalI = abs(endPhi[0])*abs(endPhi[0]) + abs(endPhi[1])*abs(endPhi[1])
        totalI = absI0_backgr[100]*absI0_backgr[100] + absIg_backgr[100]*absIg_backgr[100]

        # If more than 99% of initial intensity is lost
        # consider it reached max penetration depth
        if totalI < 0.01:
            print "max depth", depth + tiltCorrect, "Xi.", totalI
            reachedMaxDepth = True

            # Save final, tilt-corrected, depth values.
            Depth = (depth + tiltCorrect) * np.sin(-np.radians(tiltS))

        # Reinitialise the initial condition with calulated one
        phi0 = np.array([phi_backgr[100,0], phi_backgr[100,1]])

        # Move one step lower
        depth += dZ


    print "Background is at",  SumI0_backgr, SumIg_backgr
    #return [SumI0_backgr, SumIg_backgr, Depth]
    return [absI0_backgr.sum(), absIg_backgr.sum()]



def thetaB(g, a, c, V):
    from crystalography import gR_hex
    g_Matrix =  sp.Matrix([g[0], g[1], g[2]])
    d = ((g_Matrix.T * gR_hex(a, c) * g_Matrix)**0.5)[0,0]
    #constants here
    keV_J = 1.60218e-16
    Plank_c = 6.62607004e-34  # J*s
    m_e = 9.10938356e-31 # kg

    wavelength = Plank_c*(10**9)/((2.* m_e * V*keV_J)**0.5) #nm
    thetaB = wavelength/(2.*d)
    return thetaB
