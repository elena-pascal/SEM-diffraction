
import sympy as sp
from sympy.utilities.lambdify import lambdify



def defineBetaScrew(a, c, b, g, thetaB_fun, V, tiltS, rot_c):
    '''
    Return a numerical beta function dependent on x, y, z
    '''
    from sympy.physics.vector import ReferenceFrame
    from calculateBeta2 import ScrewBeta

    LF_E = ReferenceFrame('LF_E')
    thetaB = thetaB_fun(g, a, c, V)
    betaScrew_obj = ScrewBeta(a, c, b, thetaB, g, tiltS, rot_c, LF_E)
    Beta_val = betaScrew_obj.Beta_hex_RincF()

    x = sp.Symbol('x')
    y = sp.Symbol('y')
    z = sp.Symbol('z')

    return lambdify((x, y, z), Beta_val, "numpy")



def defineBetaEdge(a, c, b, irot_b, nu, g, thetaB_fun, V, tiltS, rot_c):
    '''
    Return a numerical beta function dependent on x, y, z
    '''
    from sympy.physics.vector import ReferenceFrame
    from calculateBeta2 import EdgeBeta

    LF_E = ReferenceFrame('LF_E')
    thetaB = thetaB_fun(g, a, c, V)
    betaEdge_obj = EdgeBeta(a, c, b, irot_b, nu, thetaB, g, tiltS, rot_c, LF_E)
    Beta_val = betaEdge_obj.Beta_hex_RincF()

    x = sp.Symbol('x')
    y = sp.Symbol('y')
    z = sp.Symbol('z')

    return lambdify((x, y, z), Beta_val, "numpy")



def defineBetaZero():
    '''
    Return a numerical beta function dependent on x, y, z
    '''

    x = sp.Symbol('x')
    y = sp.Symbol('y')
    z = sp.Symbol('z')

    betaPerfect = 0*x + 0*y + 0*z

    return lambdify((x, y, z), betaPerfect, "numpy")
