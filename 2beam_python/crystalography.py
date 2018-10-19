# all useful geoemtry manipulation matrixes go here
import sympy as sp
import mpmath as mp
import numpy as np
a, c = sp.symbols('a c')

########### Cubic system
def g_cubic(a):
    ''' Real space metric tensor for cubic system
    defined by lattice parameter a'''
    return sp.Matrix([[a*a, 0., 0.], [0., a*a, 0.], [0., 0., a*a]])

def gR_cubic(a):
    ''' Reciprocal space metric tensor for cubic system
    defined by lattice parameter a'''
    return sp.Matrix([[1./(a*a), 0., 0.], [0., 1./(a*a), 0], [0., 0., 1./(a*a)]])

def a_cubic(a):
    '''
    The direct space structure matrix of a cubic system
    defined by lattice parameter a
    '''
    return sp.Matrix([[a, 0., 0.], [0., a, 0.], [0., 0., a]])

def b_cubic(a):
    '''
    The reciprocal space structure matrix of a cubic system
    defined by lattice parameter a
    '''
    return a_cubic(a) * gR_cubic(a)


############# Hexagonal system
def g_hex(a, c):
    '''
    Real space metric tensor for hexagonal system
    define by lattice parameters a and c
    '''
    return sp.Matrix([[a*a, -a*a/2., 0.], [-a*a/2., a*a, 0.], [0., 0., c*c]])

def gR_hex(a, c):
    '''
    Reciprocal space metric tensor for hexagonal system
    define by lattice parameters a and c
    '''
    return sp.Matrix([[4./(3.*a*a), 2./(3.*a*a), 0.], \
                      [2./(3.*a*a), 4./(3.*a*a), 0.], \
                      [0.         , 0.         , 1./(c*c)]])

def a_hex(a, c):
    '''
    The direct space structure matrix of a hexagonal system
    define by lattice parameters a and c
    '''
    return sp.Matrix([[a, -a/2., 0.], [0., sp.sqrt(3.)*a/2., 0.], [0., 0., c]])

def b_hex(a, c):
    '''
    The reciprocal space structure matrix of a hexagonal system
    define by lattice parameters a and c
    '''
    return a_hex(a, c) * gR_hex(a, c)


def angle_between_hex(v1, v2, a, c):
    '''
    Calculates angle between two vectors v1 and v2 defined in real space
    in a hexagonal crystal frame
    '''
    # formula from Structure of Materials by M. De Graef & M. E. McHenry pg. 85
    v1v2 = sp.Matrix([ [v1[0], v1[1], v1[2]], \
                       [v2[0], v2[1], v2[2]] ])

    v1v2Dotv1v2 = sp.Matrix( v1v2 * g_hex(a, c) * v1v2.T)

    theta = sp.acos( v1v2Dotv1v2[0, 1] / (np.sqrt( np.double(v1v2Dotv1v2[0, 0]*v1v2Dotv1v2[1, 1]) ) ))

    return mp.degrees(theta)
