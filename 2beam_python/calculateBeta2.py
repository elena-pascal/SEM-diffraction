# Produce strain fields expressions for an edge dislocation in infinite medium and
# in a half infinite sphere z<0 taking into account surface relaxation
# and in a slab of material
# The x, y, z is the crystal cartesian coordinate system
# defined such that e_x//a , e_z//c* , e_y = e_z*e_x
#
from sympy.utilities.codegen import codegen
from sympy.utilities.lambdify import lambdify
from coordinates import sampleFrame, crystalFrame, displacementFrame, incidentBeamFrame

from mpmath import radians
from sympy import Eq, Matrix, atan, log, symbols, Symbol, lambdify, diff, sqrt, pi, N
from sympy.physics.vector import ReferenceFrame, express

from crystalography import b_hex, b_cubic



betaH = Symbol('betaH')
a, c = symbols('a c')
thetaB, nu, b, depth = symbols('thetaB nu b depth')
x, y, z = symbols('x y z')
rot_b = Symbol('rot_b')
rot_c = Symbol('rot_c')
tiltS = Symbol('tiltS')
tiltD = Symbol('tiltD')


DetF = ReferenceFrame('DetF')



################################################################################
#           Edge Dislocation Strain
################################################################################


class EdgeBeta(object):
    '''
    calculate beta field for an edge dislocation for a hexagonal crystal
    displacement fields are adapted from Yoffee, 1961
    tilt axis is crystal's [100] and specimen normal is along crystalographic [001]
    g is defined by the h, k, l
    return a function in right handed Cartesian system Lab frame x, y, z
    '''

    #LF = ReferenceFrame('LF')
    SF = ReferenceFrame('SF')
    CF = ReferenceFrame('CF')
    DF = ReferenceFrame('DF')
    RincF = ReferenceFrame('RincF')

    def __init__(self, a, c, b, rot_b, nu, thetaB, g, tiltS, rot_c, LF):
        self.a = a
        self.c = c
        self.b = b                                      # b = a [100] in crystal frame
        self.rot_b = rot_b                              # rotation of burger vector around dislocation line
        self.nu = nu
        self.thetaB = thetaB
        self.g = Matrix([g[0], g[1], g[2]])             # g vector h, k, l coordinates
        self.tiltS = tiltS                              # tilt of sample around [100]
        self.rot_c = rot_c                              # rotation of crystal from [100] being along the tilt axis
        self.description = 'calculate edge beta for a given configuration'
        # define required relationships between reference frames
        self.LF = LF
        # sample is tilted by tiltS wrt to the lab frame
        self.SF = sampleFrame(radians(tiltS), self.LF)
        # crystal is rotated by rot_c wrt to sample frame
        self.CF = crystalFrame(radians(rot_c), 0., self.SF)
        # displacement frame is rotated by rot_b wrt crystal frame
        self.DF = displacementFrame(radians(rot_b), self.CF)
        # the incident beam frame is rotated 180 around x wrt to the lab frame
        self.RincF = incidentBeamFrame(self.LF)
        #use the inverse of pi as a constant for computational optimisation
        self.ipi = 1./pi


    def r(self):
        '''
        This is the position vector magnitude defined in the displacement field frame
        '''
        #return ((self.DF[0]**2)+(self.DF[1]**2)+(self.DF[2]**2))**0.5
        return ((x*x)+(y*y)+(z*z))**0.5


    #define all infinte displacements
    def uInf_x(self, x, y):
        return (self.b*self.ipi/2.) * ( atan(y/x) + ( x*y/ \
            ( 2.*(1.-self.nu)*(x*x+(y*y)) ) ) )

    def uInf_y(self, x, y):
        return (-self.b*self.ipi/2.) * ( (1.-(2.*self.nu))*log(x*x+ \
            (y*y))/(4.*(1.-self.nu)) + ((x*x-(y*y)) \
            /(4.*(1.-self.nu)*(x*x+(y*y))) ))


    #use Yoffee surface relaxation formulas
    def uSurf_x(self, x, y, z):
        return (self.nu*self.b*self.ipi/(4.*(1.-self.nu))) * ( 2.*x*y*z \
            /(self.r()*((self.r() - z)**2)) + ( (1.-(2.*self.nu)) \
            *x*y/((self.r() - z)**2) ) )

    def uSurf_y(self, x, y, z):
        return (self.nu*self.b*self.ipi/(4.*(1.-self.nu))) * ( (1.-(2.*self.nu))* \
            log(self.r() - z) - ((3.-(2.*self.nu))*z/(self.r() - z)) \
            + ((3.-(2.*self.nu))*y*y/((self.r() - z)**2)) \
            - (2.*y*y/(self.r()*(self.r() - z))) )

    def uSurf_z(self, x, y, z):
         return (self.nu*self.b*y*self.ipi/(2.*(1.-self.nu))) \
            * ( (1./self.r()) + ((1.-(2.*self.nu))/(self.r() - z)) )

        #calculate total displacements
        #ECCI: infinite dislocation plus top surface relaxation
    def uTotal_x_ECCI(self, x, y, z):
        return self.uInf_x(x, y) + self.uSurf_x(x, y, z)

    def uTotal_y_ECCI(self, x, y, z):
        return self.uInf_y(x, y) + self.uSurf_y(x, y, z)

    def uTotal_z_ECCI(self, x, y, z):
        return  self.uSurf_z(x, y, z)

    # finally define a matrix with 3 entries for the x, y, z displacement components
    def u_M_RH(self, x, y, z):
        '''store the three dimensions of the displacement field in a matrix
        '''
        uM = Matrix([self.uTotal_x_ECCI(x, y, z), self.uTotal_y_ECCI(x, y, z), self.uTotal_z_ECCI(x, y, z)])
        return uM

    # u_field is defined in the awkward left handed z down coordinate frame
    # we transform this frame to a RH cartesian frame by swaping x and y
    # the u_field was defined by (rot_from_x, bx), now defined by (rot_from_y, by)
    #def u_M_RH(self, x, y, z):
    #    uM_RH = self.u_M_LH(x, y, z).subs([(x, y), (y, x)], simultaneous=True)
    #    return uM_RH

    # let's actually define it in the dislocation/displacement coordinate frame
    def u_DF(self):
        # replace x, y, z, with the reference frames base scalars R[0], R[1], R[2]
        # the vector u(ux, uy, uz) = ux*unit_x + uy*unit_y + uz*unit_z = ux*R.x + uy*R.y + uz*R.z
        uM = self.u_M_RH(x, y, z).subs([(x,self.DF[0]), (y, self.DF[1]), (z, self.DF[2])])
        uM_DF = uM[0]*self.DF.x + uM[1]*self.DF.y + uM[2]*self.DF.z
        return uM_DF

    def rinc_DF(self):
        '''
        the incident beam is along the lab cartesian direction (0,0,-1)
        return in vector form
        '''
        rinc_LF = -1. * self.LF.z
        rinc_DF = express(rinc_LF, self.DF)
        return rinc_DF

    # b in sample frame ----------
    def b_SF(self):
        # b = a * [100] in displacement frame
        #bDF = self.b * self.DF.y  # vector explicitely defined
        bDF = self.b * self.DF.x
        bSF = express(bDF, self.SF)
        x_b = bSF.dot(self.SF.x)
        y_b = bSF.dot(self.SF.y)
        z_b = bSF.dot(self.SF.z)
        return [x_b, y_b, z_b]

    # b in displacement frame ----------
    def b_DF(self):
        # b = a * [010] in displacement frame
        bDF = self.b * self.DF.x  # vector explicitely defined
        x_b = bDF.dot(self.DF.x)
        y_b = bDF.dot(self.DF.y)
        z_b = bDF.dot(self.DF.z)
        return [x_b, y_b, z_b]

    ################ No crystalography until this point #####################

    ################ cubic crystal ##########################################

    # g in dislocation frame -----------
    def gcubic_DF(self):
        #gDF = TM_ds * b_hex(self.a, self.c) * self.g
        gC = b_cubic(self.a) * self.g # g in real space Cartesian coordinates
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gDF = express(gCF, self.DF)
        return gDF

    # g in sample frame -----------
    def gcubic_SF(self):
        gC = b_cubic(self.a) * self.g #  g in real space Cartesian coordinates
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gSF = express(gCF, self.SF)
        x_g = gSF.dot(self.SF.x)
        y_g = gSF.dot(self.SF.y)
        z_g = gSF.dot(self.SF.z)
        return [x_g, y_g, z_g]

    #===============================================================================
    #                    Beta
    #===============================================================================
    def Beta_cubic_DF(self):
        '''
        b_curvature = (r_inc^d)_i * d(u^d)_j/dx_i * (T^ds)_jl * B_lm * g_m
        b_displacement = grad(u_i * a_hex_il * gR_hex_lj * gR_j)*r_g(x, y, z)
        a_hex_il * gR_hex_lj = b_hex_ij
        b_tunstall = b_curvature + thetaB*b_displacemnt
        Return result in displacement frame
        '''



        udotg = self.u_DF().dot(self.gcubic_DF())

        b_curv =  diff(udotg, self.DF[0]) * (self.rinc_DF().dot(self.DF.x))  + \
                  diff(udotg, self.DF[1]) * (self.rinc_DF().dot(self.DF.y))  + \
                  diff(udotg, self.DF[2]) * (self.rinc_DF().dot(self.DF.z))

        rg = self.gcubic_DF().normalize()

        b_displ =  diff(udotg, x) * (rg.dot(self.DF.x))  + \
                   diff(udotg, y) * (rg.dot(self.DF.y))  + \
                   diff(udotg, z) * (rg.dot(self.DF.z))

        betaT_DF = b_curv + self.thetaB * b_displ

        #return express(betaT_DF, self.SF, variables = True)
        return betaT_DF

    def Beta_cubic_DF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting'''
        beta = self.Beta_cubic_DF().subs([(self.DF[0], x), (self.DF[1], y), (self.DF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), beta, "numpy")
        return np_function(xN, yN, zN)

    def Beta_cubic_LF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in lab frame'''
        betaLF = express(self.Beta_cubic_DF(), self.LF, variables = True)
        betaLF_n = betaLF.subs([(self.LF[0], x), (self.LF[1], y), (self.LF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaLF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_cubic_SF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in sample frame'''
        betaSF = express(self.Beta_cubic_DF(), self.SF, variables = True)
        betaSF_n = betaSF.subs([(self.SF[0], x), (self.SF[1], y), (self.SF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaSF_n, "numpy")
        return np_function(xN, yN, zN)





    ################ hexagonal crystal ######################################
    # g in dislocation frame -----------
    def ghex_DF(self):
        gC = b_hex(self.a, self.c) * self.g # g in real space Cartesian coordinates
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gDF = express(gCF, self.DF)
        return gDF

    # g in sample frame -----------
    def ghex_SF(self):
        gC = b_hex(self.a, self.c) * self.g #  g in real space Cartesian coordinates
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gSF = express(gCF, self.SF)
        x_g = gSF.dot(self.SF.x)
        y_g = gSF.dot(self.SF.y)
        z_g = gSF.dot(self.SF.z)
        return [x_g, y_g, z_g]


    def Beta_hex_DF(self):
        '''
        b_curvature = (r_inc^d)_i * d(u^d)_j/dx_i * (T^ds)_jl * B_lm * g_m
        b_displacement = grad(u_i * a_hex_il * gR_hex_lj * gR_j)*r_g(x, y, z)
        a_hex_il * gR_hex_lj = b_hex_ij
        b_tunstall = b_curvature + thetaB*b_displacemnt
        Return result in displacement frame
        '''

        udotg = self.u_DF().dot(self.ghex_DF())
        #udotg_LB = express(udotg, LB, variables=True)

        # directional derivative in displacement frame
        # find better way to get projections on basis vectors?
        b_curv = ( diff(udotg, self.DF[0]) * (self.rinc_DF().dot(self.DF.x)) ) + \
                 ( diff(udotg, self.DF[1]) * (self.rinc_DF().dot(self.DF.y)) ) + \
                 ( diff(udotg, self.DF[2]) * (self.rinc_DF().dot(self.DF.z)) )

        rg = self.ghex_DF().normalize()

        b_displ = ( diff(udotg, x) * (rg.dot(self.DF.x)) ) + \
                  ( diff(udotg, y) * (rg.dot(self.DF.y)) ) + \
                  ( diff(udotg, z) * (rg.dot(self.DF.z)) )

        betaT_DF = b_curv + (self.thetaB * b_displ)
        #return express(betaT_DF, self.SF, variables = True)
        return betaT_DF

    def Beta_hex_DF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting'''
        beta = self.Beta_hex_DF().subs([(self.DF[0], x), (self.DF[1], y), (self.DF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), beta, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_LF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in lab frame'''
        betaLF = express(self.Beta_hex_DF(), self.LF, variables = True)
        betaLF_n = betaLF.subs([(self.LF[0], x), (self.LF[1], y), (self.LF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaLF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_SF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in sample frame'''
        betaSF = express(self.Beta_hex_DF(), self.SF, variables = True)
        betaSF_n = betaSF.subs([(self.SF[0], x), (self.SF[1], y), (self.SF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaSF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_RincF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in incident beam frame'''
        betaRincF = express(self.Beta_hex_DF(), self.RincF, variables = True)
        betaRincF_n = betaRincF.subs([(self.RincF[0], x), (self.RincF[1], y), (self.RincF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaRincF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_RincF(self):
        ''' numerical function of above for plotting
            Return result in incident beam frame'''
        betaRincF = express(self.Beta_hex_DF(), self.RincF, variables = True)
        betaRincF_n = betaRincF.subs([(self.RincF[0], x), (self.RincF[1], y), (self.RincF[2], z)], simultaneous=True)
        return betaRincF_n
#===============================================================================

#======================== write f95 module =====================================
    def writeModule(self, fname):
        expr = Eq(betaH, self.Beta_hex_RincF_n(x, y, z))
        fct = codegen(("beta_ECCI", expr), "f95", fname, to_files = True,\
                        argument_sequence =(x, y, z, betaH),\
                        project ='Edgebeta_generation')




################################################################################
#           Screw Dislocation Strain
################################################################################

class ScrewBeta(object):
    '''
    calculate beta field for an screw dislocation for different crystal configs
    displacement fields are adapted from Yoffee, 1961
    return a formula depending on variables x, y, z
    '''

    #LF = ReferenceFrame('LF')
    SF = ReferenceFrame('SF')
    CF = ReferenceFrame('CF')
    DF = ReferenceFrame('DF')
    RincF = ReferenceFrame('RincF')


    def __init__(self, a, c, b, thetaB, g, tiltS, rot_c, LF):
        self.a = a
        self.c = c
        self.b = b
        self.thetaB = thetaB
        self.g = Matrix([g[0], g[1], g[2]])             # g vector h, k, l coordinates
        self.tiltS = tiltS                              # tilt of sample around [100]
        self.rot_c = rot_c                              # rotation of crystal from [100] being along the tilt axis
        self.description = 'calculate screw beta for a given configuration'
        # define required relationships between reference frames

        # define required relationships between reference frames
        self.LF = LF
        # sample is tilted by tiltS wrt to the lab frame
        self.SF = sampleFrame(radians(tiltS), self.LF)
        # crystal is rotated by rot_c wrt to sample frame
        self.CF = crystalFrame(radians(rot_c), 0., self.SF)
        # displacement frame is rotated by rot_b wrt crystal frame
        self.DF = self.CF# displacementFrame is the same as crystalFrame for a screw dislocation
        # the incident beam frame is rotated 180 around x wrt to the lab frame
        self.RincF = incidentBeamFrame(self.LF)


    def r(self):
        '''
        This is the position vector magnitude defined in the displacement field frame
        '''
        #return ((DF[0]**2)+(DF[1]**2)+(DF[2]**2))**0.5
        return ((x*x)+(y*y)+(z*z))**0.5

    #define all infinte displacements
    def uInf_z(self, x, y):
        return (-self.b/(2.*pi)) *  atan(y/x)


    #use Yoffee surface relaxation formulas
    def uSurf_x(self,  z):
        return (self.b/(2.*pi)) * y/( self.r() - z )

    def uSurf_y(self,  z):
        return (-self.b/(2.*pi)) * x/( self.r() - z )



        #calculate total displacements
        #ECCI: infinite dislocation plus top surface relaxation
    def uTotal_x_ECCI(self, x, y, z):
        return self.uSurf_x(z)

    def uTotal_y_ECCI(self, x, y, z):
        return self.uSurf_y(z)

    def uTotal_z_ECCI(self, x, y):
        return self.uInf_z(x, y)

    # finally define a matrix with 3 entries for the x, y, z displacement components
    def u_M_RH(self, x, y, z):
        '''store the three dimensions of the displacement field in a matrix
        '''
        uM = Matrix([self.uTotal_x_ECCI(x, y, z), self.uTotal_y_ECCI(x, y, z), self.uTotal_z_ECCI(x, y)])
        #uM = Matrix([0, 0, self.uTotal_z_ECCI(x, y)])
        return uM

    # u_field is defined in the awkward left handed z down coordinate frame
    # we transform this frame to a RH cartesian frame by swaping x and y
    # the u_field was defined by (rot_from_x, bx), now defined by (rot_from_y, by)
    #def u_M_RH(self, x, y, z):
    #    uM_RH = self.u_M_LH(x, y, z).subs([(x, y), (y, x)], simultaneous=True)
    #    return uM_RH

    # let's actually define it in the dislocation/displacement coordinate frame
    def u_DF(self):
        # replace x, y, z, with the reference frames base scalars R[0], R[1], R[2]
        # the vector u(ux, uy, uz) = ux*unit_x + uy*unit_y + uz*unit_z = ux*R.x + uy*R.y + uz*R.z
        uM = self.u_M_RH(x, y, z).subs([(x,self.DF[0]), (y, self.DF[1]), (z, self.DF[2])])
        uM_DF = uM[0]*self.DF.x + uM[1]*self.DF.y + uM[2]*self.DF.z
        return uM_DF

    def rinc_DF(self):
        '''
        the incident beam is along the lab cartesian direction (0,0,-1)
        return in vector form
        '''
        rinc_LF = -1. * self.LF.z
        rinc_DF = express(rinc_LF, self.DF)
        return rinc_DF

    # g in dislocation frame -----------
    def g_DF(self):
        #gDF = TM_ds * b_hex(self.a, self.c) * self.g
        gC = b_hex(self.a, self.c) * self.g # still matrix form
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gDF = express(gCF, self.DF)
        return gDF

    # g in sample frame -----------
    def g_SF(self):
        gC = b_hex(self.a, self.c) * self.g #  g in real space Cartesian coordinates
        gCF = gC[0]* self.CF.x + gC[1]*self.CF.y + gC[2]*self.CF.z # vector explicitely defined
        gSF = express(gCF, self.SF)
        x_g = gSF.dot(self.SF.x)
        y_g = gSF.dot(self.SF.y)
        z_g = gSF.dot(self.SF.z)
        return [x_g, y_g, z_g]

    # b in sample frame ----------
    def b_SF(self):
        bDF = self.b * self.DF.z
        bSF = express(bDF, self.SF)
        x_b = bSF.dot(self.SF.x)
        y_b = bSF.dot(self.SF.y)
        z_b = bSF.dot(self.SF.z)
        return [x_b, y_b, z_b]
    #===============================================================================
    #                    beta
    #===============================================================================
    def Beta_hex_DF(self):
        '''
        b_curvature = (r_inc^d)_i * d(u^d)_j/dx_i * (T^ds)_jl * B_lm * g_m
        b_displacement = grad(u_i * a_hex_il * gR_hex_lj * gR_j)*r_g(x, y, z)
        a_hex_il * gR_hex_lj = b_hex_ij
        b_tunstall = b_curvature + thetaB*b_displacemnt
        Return result in displacement frame
        '''
        udotg = self.u_DF().dot(self.g_DF())
        #udotg_LB = express(udotg, LB, variables=True)

        # directional derivative in displacement frame
        # find better way to get projections on basis vectors?
        b_curv = diff(udotg, self.DF[0]) * (self.rinc_DF().dot(self.DF.x)) + \
                 diff(udotg, self.DF[1]) * (self.rinc_DF().dot(self.DF.y)) + \
                 diff(udotg, self.DF[2]) * (self.rinc_DF().dot(self.DF.z))

        #rgC = gC_L / sqrt(gC_L[0]*gC_L[0] + (gC_L[1]*gC_L[1]) + (gC_L[2]*gC_L[2]))
        rg = self.g_DF().normalize()

        b_displ = diff(udotg, x) * (rg.dot(self.DF.x)) + \
            diff(udotg, y) * (rg.dot(self.DF.y)) + \
            diff(udotg, z) * (rg.dot(self.DF.z))

        betaT_DF = b_curv + self.thetaB * b_displ
        #return express(betaT_DF, self.SF, variables = True)
        return betaT_DF

    def Beta_hex_DF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting'''
        beta = self.Beta_hex_DF().subs([(self.DF[0], x), (self.DF[1], y), (self.DF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), beta, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_LF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in lab frame'''
        betaLF = express(self.Beta_hex_DF(), self.LF, variables = True)
        betaLF_n = betaLF.subs([(self.LF[0], x), (self.LF[1], y), (self.LF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaLF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_SF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in sample frame'''
        betaSF = express(self.Beta_hex_DF(), self.SF, variables = True)
        betaSF_n = betaSF.subs([(self.SF[0], x), (self.SF[1], y), (self.SF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaSF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_RincF_n(self, xN, yN, zN):
        ''' numerical function of above for plotting
            Return result in incident beam frame'''
        betaRincF = express(self.Beta_hex_DF(), self.RincF, variables = True)
        betaRincF_n = betaRincF.subs([(self.RincF[0], x), (self.RincF[1], y), (self.RincF[2], z)], simultaneous=True)
        np_function = lambdify((x, y, z), betaRincF_n, "numpy")
        return np_function(xN, yN, zN)

    def Beta_hex_RincF(self):
        ''' numerical function of above for plotting
            Return result in incident beam frame'''
        betaRincF = express(self.Beta_hex_DF(), self.RincF, variables = True)
        betaRincF_n = betaRincF.subs([(self.RincF[0], x), (self.RincF[1], y), (self.RincF[2], z)], simultaneous=True)
        return betaRincF_n

#======================== write f95 module =====================================
    def writeModule(self, fname):
        expr = Eq(betaH, self.Beta_hex_RincF_n(x, y, z))
        fct = codegen(("beta_ECCI", expr), "f95", fname, to_files = True,\
                        argument_sequence =(x, y, z, betaH),\
                        project ='Screwbeta_generation')
