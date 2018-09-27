from mpmath import radians
from sympy.physics.vector import ReferenceFrame, express
from sympy import Symbol, symbols, Eq
from coordinates import sampleFrame, crystalFrame, displacementFrame, detectorFrame
from sympy.utilities.codegen import codegen
from calculateBeta2 import EdgeBeta
from calculateBeta2 import ScrewBeta

x, y, z = symbols('x y z')
LFe = ReferenceFrame('LFe')
LFs = ReferenceFrame('LFs')
LFm = ReferenceFrame('LFm')



### diffraction conditions
g1 = (7.,5.,-3.)
g2 = (-5.,-7.,-3.)

# material parameters
a = 0.32   #GaN
c = 0.52

# Edge Burgers vectors--------------------------
be = a
# dictionary containing Burgers vectors possible rotation angles in degrees
be_vectors = {}
for i in range(0, 6):
    be_vectors['rot_b%01d'%i] = i * 60 # in degrees

#rot_b1 = 0    # [100] in Cartesian coordinates
#rot_b2 = 60   # [1-10]
#rot_b3 = 120  # [0-10]
#rot_b4 = 180  # [-100]
#rot_b5 = 240  # [-110]
#rot_b6 = 300  # [010]

# Screw Burgers vectors--------------------------
bs = c


# Mixed Burgers vectors--------------------------
bm = (a*a + c*c)**0.5



tiltS = -49.6                          # tilt of sample around [100]

tiltD = 20                             # detector normal tilt from horizontal

rot_c1 = 58.35                         # crystal rotation
rot_c2 = 61.65

nu = 0.3
thetaB = 0.005



##################################################
#   6 Edge dislocation strain modules per g vector
##################################################

#def __init__(self, a, c, b, rot_b, nu, thetaB, g, tiltS, rot_c, LF)

for key, irot_b in sorted(be_vectors.items()):

    betaEdge_hex1 = EdgeBeta(a, c, be, irot_b, nu, thetaB, g1, tiltS, rot_c1, LFe)
    betaEdge_hex1.writeModule('beta_modules/betaEdgeHex_g1_b%s' %irot_b)

    betaEdge_hex2 = EdgeBeta(a, c, be, irot_b, nu, thetaB, g2, tiltS, rot_c2, LFe)
    betaEdge_hex2.writeModule('beta_modules/betaEdgeHex_g2_b%s' %irot_b)

    print "Wrote Edge Beta modules g1, g2 for rot_b = ", irot_b

print

##################################################
#  2 Screw dislocation strain modules per g vector
##################################################

#def __init__(self, a, c, b, thetaB, g, tiltS, rot_c, LF):

for sign in [-1., 1.,]:
    betaScrew_hex1 = ScrewBeta(a, c, sign*bs, thetaB, g1, tiltS, rot_c1, LFs)
    betaScrew_hex1.writeModule('beta_modules/betaScrewHex_g1_b%s' %sign)

    betaScrew_hex2 = ScrewBeta(a, c, sign*bs, thetaB, g2, tiltS, rot_c2, LFs)
    betaScrew_hex2.writeModule('beta_modules/betaScrewHex_g2_b%s' %sign)

    print "Wrote Screw Beta modules g1, g2 for bs = ", sign, "c"

print

####################################################
#  12 Mixed dislocation strain modules per g vector
####################################################

betaH = Symbol('betaH')
for sign in [-1., 1.]:
    for key, irot_b in sorted(be_vectors.items()):

        betaEdge_hex1 = EdgeBeta(a, c, be, irot_b, nu, thetaB, g1, tiltS, rot_c1, LFm)
        strainEdge = betaEdge_hex1.Beta_hex_RincF_n( x, y, z)

        betaScrew_hex1 = ScrewBeta(a, c, sign*bs, thetaB, g1, tiltS, rot_c1, LFm)
        strainScrew = betaScrew_hex1.Beta_hex_RincF_n( x, y, z)


        expr = Eq(betaH, strainEdge + strainScrew)

        fct = codegen(("beta_ECCI", expr), "f95", 'beta_modules/betaMixedHex_g1_bs%01d_be%s' %(sign, irot_b), to_files = True,\
                argument_sequence =(x, y, z, betaH),\
                project ='mixedbeta_generation')



        betaEdge_hex2 = EdgeBeta(a, c, be, irot_b, nu, thetaB, g2, tiltS, rot_c2, LFm)
        strainEdge = betaEdge_hex1.Beta_hex_RincF_n( x, y, z)

        betaScrew_hex2 = ScrewBeta(a, c, sign*bs, thetaB, g2, tiltS, rot_c2, LFm)
        strainScrew = betaScrew_hex1.Beta_hex_RincF_n( x, y, z)

        expr = Eq(betaH, strainEdge + strainScrew)
        fct = codegen(("beta_ECCI", expr), "f95", 'beta_modules/betaMixedHex_g2_bs%01d_be%s' %(sign, irot_b), to_files = True,\
                argument_sequence =(x, y, z, betaH),\
                project ='mixedbeta_generation')

        print "Wrote Mixed Beta modules g1, g2 for bs = ", sign, "c, be = ", irot_b
