# Defines relationship between the different coordinate systems.
# The vector module of sympy physics assumes RH coordinate systems
# therefore, all reference frames defined below will be right handed
# All rotations are defined to be anticlockwise
from sympy.physics.vector import ReferenceFrame, express
from mpmath import radians
import sympy as sp
gx, gy, gz = sp.symbols('gx gy gz')
tiltS = sp.Symbol('tiltS') # sample tilt from horizontal
tiltD = sp.Symbol('tiltD') # angle between detector normal and horizontal
c_rot = sp.Symbol('c_rot')
b_rot = sp.Symbol('b_rot')
z_tilt = sp.Symbol('z_tilt')

'''
The static reference frame will be the lab frame
'''
LF = ReferenceFrame('LF')

SF = ReferenceFrame('SF')
detF = ReferenceFrame('detF')
CF = ReferenceFrame('CF')
RincF = ReferenceFrame('RincF')
DF = ReferenceFrame('DF')

def sampleFrame(tiltS, LF):
    '''
    The sample reference frame is tilted wrt to lab frame
    around the x axis assuming LB.x//SF.x by an angle tiltS
    LF is the lab frame and must be of type ReferenceFrame
    '''
    SF.orient(LF, 'Axis', [tiltS, LF.x])
    return SF

def detectorFrame(tiltD, LF):
    '''
    The sample reference frame is tilted in the lab frame
    such that the angle between lab.y and the detector normal is tiltD
    to go from the lab frame to the detector frame:
    1) rotate lab frame around lab.z by angle pi
    2) rotate new frame around new x by angle tiltD-pi/2
    '''
    detF.orient(LF, 'Body', [radians(180), tiltD - radians(90), 0], 'ZXZ')
    return detF


def crystalFrame(c_rot, z_tilt, SF):
    '''
    The crystal Cartesian reference frame is rotated wrt to the sample frame
     1) such that it meets the crystal definition of x being the tilt Axis
     2) then rotated about z to account for the crystal rotations
     3) tilted around x to eventually relax CF.z//SF.z (?)

    SF is the sample frame and must be of type ReferenceFrame
    '''
    CF.orient(SF, 'Axis', [c_rot, SF.z])

    #CF_tilted = ReferenceFrame('CF_tilted')
    #CF_tilted.orient(CF, 'Axis', [z_tilt, CF.z])
    return CF



def displacementFrame(b_rot, CF):
    '''
     The edge dislocation displacement field left handed cartesian system (z down)
     is defined by the Burgers' vector b_x
     When we transform this system to a right handed system (z down) by swaping x and y
     the defining Burgers' vector will be along y such that DF.y // CF.x
     The displacement field reference frame (DF) is defined in the crystal frame by
          rotation around CF.z by rot_b
    CF is the crystal frame and must be of type ReferenceFrame

     Changed the old version above to Yoffee's system. Displcement frame
     is now RH with z up denfined by a rot_b rotation around z
    '''

    DF.orient(CF, 'Axis', [b_rot,CF.z])
    return DF



def incidentBeamFrame(LF):
    '''
    Eventually, the strain profile is to be saved in r_inc frame
    which is just the lab frame roated by 180 around x
    '''

    RincF.orient(SF, 'Axis', [radians(180),SF.x])
    return RincF


# The incident beam vector is defined in the lab frame system by
incBeam_LF = -1. * LF.z
