import numpy as np
from numpy import pi

def zfunc(y, t, frac0, fracg, w, beta):
    '''
    Where:
     y = initial values
     frac0 = RPAR(2) = Xsi0'/Xsig'
     fracg = RPAR(3) = Xsig/Xsig'

     Xsi0'= mean/background absorbtion acting equally on both Bloch waves
     Xsig = extinction distancefor beam g
     Xsig'= anomalous absorbtion distance

    '''
    y1, y2 = y
    fg_f0 = fracg/frac0

    #YDOT(1) = -PI*Y(1)*RPAR(3)/RPAR(2) + PI*Y(2)*(im-RPAR(3))
    #YDOT(2) = PI*Y(1)*(im-RPAR(3)) + PI*Y(2)*(2.0_dp*im*RPAR(1)-(RPAR(3)/RPAR(2))+2.0_dp*im*beta)
    y1dot = -pi*y1*fg_f0 + pi*y2*(1j-fracg)
    y2dot = pi*y1*(1j-fracg) + pi*y2*(2.*1j*w - fg_f0 + 2.*1j*beta)

    return [y1dot, y2dot]

# The Jacobian is
def zjac(y, t, frac0, fracg, w, beta):
    y1, y2 = y
    fg_f0 = fracg/frac0

    j11 = -pi*fg_f0
    j12 = pi*(1j-fracg)
    j21 = pi*(1j-fracg)
    j22 = pi*(2.*1j*w - fg_f0 + 2.*1j*beta)

    jac = np.array([[j11, j12], [j21, j22]])
    return jac
