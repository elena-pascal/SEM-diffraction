# Copyright (c) Elena Pascal Physics, University of Strathclyde.
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3.




import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sympy as sp
import numpy.ma as ma
import mpmath as mp
import mayavi.mlab as mlab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator

from sympy import Matrix
from sympy.physics.vector import ReferenceFrame, express


nu, b, thetaB = sp.symbols('nu b thetaB')
a, c = sp.symbols('a c')
x, y, z = sp.symbols('x y z')
rotb = sp.Symbol('rotb')
rotc = sp.Symbol('rotc')
tiltS = sp.Symbol('tiltS')
tiltD = sp.Symbol('tiltD')

# Why do I need different frames? aka the second time I call the same frame I get no relationship
# edge frames
LF_E = ReferenceFrame('LF_E')

# screw frames
LF_S = ReferenceFrame('LF_S')

from crystalography import b_hex, b_cubic, gR_hex, angle_between_hex
from calculateBeta2 import EdgeBeta, ScrewBeta

from fileTools import readInput, writeOutput, fileName

# Read parameters from file.
my_file = "plotStrain.in"
indata = readInput(my_file)

a = indata['a'] #nm
c = indata['c'] #nm

tiltS = indata['tiltS']                # tilt of sample around [100]x
tiltD = indata['tiltD']                # detector normal tilt from horizontal

g = indata['g']

rot_c = indata['rot_c']                 # crystal rotation

nuGaN = indata['nu']

# Calculate angle between g vectors
#g1R = b_hex(a, c) * Matrix([g1[0], g1[1], g1[2]])
#g2R = b_hex(a, c) * Matrix([g2[0], g2[1], g2[2]])
#print angle_between_hex(g1R, g2R, a, c)




# Calculate theta_B
g_Matrix =  Matrix([g[0], g[1], g[2]])
d = ((g_Matrix.T * gR_hex(a, c) * g_Matrix)**0.5)[0,0]
thetaB_GaN = 0.007/(2*d)
print "thetaB", thetaB_GaN

elif indata['type'] == 'Edge':
    my_beta = defineBetaEdge(indata['a'], indata['c'], indata['be'],\
                     indata['rot_b'], indata['nu'], thetaB, indata['V'],\
                     indata['g'], indata['tiltS'], indata['rot_c'])
else:
    my_beta = 0. # perfect crystal


if indata['type'] == 'Screw':





    #### Edge Burgers vectors ########################################

    be = a

    # dictionary containing Burgers vectors possible rotation angles in degrees
    #be_vectors = {}
    #for i in range(1, 7):
    #    be_vectors['rot_b%01d'%i] = (i-1) * 60 # in degrees

    # rotation for edge dislocaiton Burgers vector
#rot_b1 = 0    # [100] in Cartesian coordinates
#rot_b2 = 60   # [1-10]
#rot_b3 = 120  # [0-10]
#rot_b4 = 180  # [-100]
#rot_b5 = 240  # [-110]
#rot_b6 = 300  # [010]
irot_b = indata['rot_b']



# here the strain field is defined
#a, c, b, rot_b, nu, thetaB, g, tiltS, rot_c, LF
betaEdge_hex = EdgeBeta(a, c, be, irot_b, nuGaN, thetaB_GaN, g, tiltS, rot_c, LF_E)


print "Formulated strain..."

# map strain on grid
xM, yM, zM = np.mgrid[-30:30:500j, -30:30:500j, -100:0.:20j] #nm
valE = betaEdge_hex.Beta_hex_SF_n(xM, yM, zM)

print
print "Mapped edge TD strain..."
#print valE


bE_inSF = betaEdge_hex.b_SF()


######### Mixed Burgers vectors ##############################################
bm = (a*a + c*c)**0.5

valM = valS + valE

######################################################################
# vector b values in sample frame
bS_inSF = betaScrew_hex.b_SF()


# vector g values in sample frame
g_inSF = betaScrew_hex.g_SF()
#g_inSF = betaEdge_hex1.ghex_SF()
#g_inSF = betaEdge_Cubic.gcubic_SF()
#g_inDF = [betaEdge_Cubic.gcubic_DF().dot(betaEdge_Cubic.DF.x),
#            betaEdge_Cubic.gcubic_DF().dot(betaEdge_Cubic.DF.y),
#            betaEdge_Cubic.gcubic_DF().dot(betaEdge_Cubic.DF.z) ]

bM_inSF = np.array(bE_inSF) + np.array(bS_inSF)
print 'g vector in SF is:', g_inSF
print 'b vector in SF is:', bM_inSF


######################################################
########### Plot here - Sample Frame #################
mlab.figure(1, size=(500, 250), bgcolor=(0, 0, 0))
mlab.clf()

# to make sure the middle colour is normalise to 0 set symmetric vmin and vmanx
cont = mlab.contour3d(xM, yM, zM, valE, contours=[-0.1, 0.1], transparent=False, colormap='RdBu', opacity = 0.9, vmin=-0.1, vmax=0.1  )
zeroline = mlab.contour3d(xM, yM, zM, valE, contours=[0.00], transparent=False, colormap='RdBu', opacity = 0.2, vmin=-0.1, vmax=0.1  )

# Transparency filter for contribution to incident beam ###############
# contribution = 10. * np.log(-zM/40.)

# Retrieve the LUT of the surf object.
#lut = cont.module_manager.scalar_lut_manager.lut.table.to_array()

# We modify the alpha channel to add a transparency gradient to the colormap :(
#lut[:, -1] = np.linspace(255, 150, 256)

# and finally we put this LUT back in the surface object. We could have
# added any 255*4 array rather than modifying an existing LUT.
#cont.module_manager.scalar_lut_manager.lut.table = lut

# reverse lut
cont.module_manager.scalar_lut_manager.reverse_lut = True
###########################################################################

#mlab.axes(cont, color=(.7, .7, .7), xlabel='', ylabel='',zlabel='',
#            z_axis_visibility=False,)

# draw some lines for b and g

bvector = mlab.plot3d(np.array([0, 100.*np.double(bM_inSF[0])]), np.array([0, 100.*np.double(bM_inSF[1])]), np.array([0., 100.*np.double(bM_inSF[2])]), color=(0.2, 0.2, 0.2), tube_radius=0.5)
gvector = mlab.plot3d(np.array([0, 1*np.double(g_inSF[0])]), np.array([0, 1*np.double(g_inSF[1])]), np.array([0., 1*np.double(g_inSF[2])]), color=(0.5, 0 ,0), tube_radius=0.5)
disline = mlab.plot3d(np.array([0, 0]), np.array([0, 0]), np.array([0., -100]), color=(1, 1 ,1), tube_radius=1)

mlab.text3d(-7, +30, -10, 'g', color=(0.5, 0 ,0), scale=5.)


#mlab.outline(cont)


# Show the colorbar bar (legend)
mlab.colorbar(title='Strain', orientation='vertical', nb_labels=3)

# Show image
#mlab.draw()
mlab.show()

# Save image
#mlab.savefig("Screw.obj")
