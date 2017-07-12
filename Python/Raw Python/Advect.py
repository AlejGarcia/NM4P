# advect - Program to solve the advection equation 
# using the various hyperbolic PDE schemes

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Select numerical parameters (time step, grid spacing, etc.).
method = input('Choose a numerical method, 1) FTCS; 2) Lax; 3) Lax-Wendroff :')
N = input('Enter number of grid points: ')
L = 1.      # System size
h = L/N     # Grid spacing
c = 1.      # Wave speed
print 'Time for wave to move one grid spacing is ', h/c 
tau = input('Enter time step: ')
coeff = -c*tau/(2.*h)    # Coefficient used by all schemes
coefflw = 2*coeff**2     # Coefficient used by L-W scheme
print 'Wave circles system in ', L/(c*tau), ' steps' 
nStep = input('Enter number of steps: ')

#* Set initial and boundary conditions.
sigma = 0.1                  # Width of the Gaussian pulse
k_wave = np.pi/sigma         # Wave number of the cosine
x = np.arange(N)*h - L/2.    # Coordinates of grid points
# Initial condition is a Gaussian-cosine pulse
a = np.empty(N)
for i in range(N) :
    a[i] = np.cos(k_wave*x[i]) * np.exp(-x[i]**2/(2*sigma**2)) 
# Use periodic boundary conditions
ip = np.arange(N) + 1  
ip[N-1] = 0          # ip = i+1 with periodic b.c.
im = np.arange(N) - 1  
im[0] = N-1          # im = i-1 with periodic b.c.

#* Initialize plotting variables.
iplot = 1           # Plot counter
nplots = 50         # Desired number of plots
aplot = np.empty((N,nplots))
tplot = np.empty(nplots)
aplot[:,0] = np.copy(a)     # Record the initial state
tplot[0] = 0                # Record the initial time (t=0)
plotStep = nStep/nplots +1  # Number of steps between plots

#* Loop over desired number of steps.
for iStep in range(nStep):  ## MAIN LOOP ##

    #* Compute new values of wave amplitude using FTCS, 
    #%  Lax or Lax-Wendroff method.
    if method == 1 :      ### FTCS method ###
        a[:] = a[:] + coeff*( a[ip] - a[im] )  
    elif  method == 2 :   ### Lax method ###
        a[:] = .5*( a[ip] + a[im] ) + coeff*( a[ip] - a[im] )   
    else:                 ### Lax-Wendroff method ###
        a[:] = ( a[:] + coeff*( a[ip] - a[im] ) + 
                coefflw*( a[ip] + a[im] -2*a[:] ) )

    #* Periodically record a(t) for plotting.
    if (iStep+1) % plotStep < 1 :        # Every plot_iter steps record 
        aplot[:,iplot] = np.copy(a)      # Record a(i) for ploting
        tplot[iplot] = tau*(iStep+1)
        iplot += 1
        print iStep, ' out of ', nStep, ' steps completed'


#* Plot the initial and final states.
plt.plot(x,aplot[:,0],'-',x,a,'--');
plt.legend(['Initial  ','Final'])
plt.xlabel('x')  
plt.ylabel('a(x,t)')
plt.show()

#* Plot the wave amplitude versus position and time
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection = '3d')
Tp, Xp = np.meshgrid(tplot[0:iplot], x)
ax.plot_surface(Tp, Xp, aplot[:,0:iplot], rstride=1, cstride=1, cmap=cm.gray)
ax.view_init(elev=30., azim=190.)
ax.set_ylabel('Position') 
ax.set_xlabel('Time')
ax.set_zlabel('Amplitude')
plt.show()

