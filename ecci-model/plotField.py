from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sympy as sp
import numpy.ma as ma
import mpmath as mp
from math import sin, cos
from matplotlib.ticker import MultipleLocator

nu, b, thetaB = sp.symbols('nu b thetaB')


from calculateBeta import EdgeBeta

a = 0.4
c = 0.4
tilt = 0
rot = 0
rot_b = 180
g_hkl = (4,0,0)
displacement = EdgeBeta(a, c, rot_b, g_hkl, tilt, rot)


X, Z = np.mgrid[-1:1:12j, 0:2:20j]
Y = -0.2
b = 1.0
nu = 0.3

#XR = X*cos(mp.radians(tilt))+ Z*sin(mp.radians(tilt))
#ZR = -X*sin(mp.radians(tilt)) + Z*cos(mp.radians(tilt))

ux = displacement.NuTotal_x_ECCI(X, Y, Z, b, nu)
uz = displacement.NuTotal_z_ECCI(X, Y, Z, b, nu)

speed = np.sqrt(ux**2 + uz**2)
uxN = ux/speed
uzN = uz/speed



plot1 = plt.figure()
#plt.axes(frameon=0)
plt.gca().invert_yaxis()
plt.quiver(X, Z, uxN, uzN,        # data
           speed,                   # colour the arrows based on this array
           cmap=cm.spectral,     # colour map
           #arrowsize=1.5,        # length of the arrows
           #pivot = 'mid',
           #arrowstyle ='->'
           #color = 'r',
           #units = 'width',
           #scale = 1/0.05
             )
plt.colorbar()                  # adds the colour bar




plt.title('Projection of the displacement field')
plt.xlabel('X[100]')
plt.ylabel('Z[001]')


#fontsize
rc('text', usetex=True)
rc('font',**{'family':'cursive','cursive':['Zapf Chancery']})
rc('font',size=18)
rc('font',family='cursive')
rc('axes',labelsize=20)


plt.gca().xaxis.set_major_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(1))


plt.savefig('edge_cubic.png')
plt.show(plot1)                 # display the plot



#plot3 =fig.set_size_inches(18.5, 10.5, forward=True) plt.figure()
#X = np.linspace(-1,1, 30)
#Z = np.linspace(-0,2, 30)
#plt.streamplot(X, Z, ux, uz, density=0.6)          # data
#               color=speed,         # array that determines the colour
#               cmap=cm.cool,        # colour map
#               linewidth=2,         # line thickness
#               arrowstyle='->')     # arrow style
#               arrowsize=1.5)       # arrow size

#plt.colorbar()                      # add colour bar on the right

#plt.title('tilt=0, rot=0, b=[010], g=[100]')
#plt.show(plot3)                     # display the plot
