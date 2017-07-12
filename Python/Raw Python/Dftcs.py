# dftcs - Program to solve the diffusion equation 
# using the Forward Time Centered Space (FTCS) scheme.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Initialize parameters (time step, grid spacing, etc.).
tau = input('Enter time step: ')
N = input('Enter the number of grid points: ')
L = 1.        # The system extends from x=-L/2 to x=L/2
h = L/(N-1)   # Grid size
kappa = 1.    # Diffusion coefficient
coeff = kappa*tau/h**2
if coeff < 0.5 :
    print 'Solution is expected to be stable'
else:
    print 'WARNING: Solution is expected to be unstable'


#* Set initial and boundary conditions.
tt = np.zeros(N)                # Initialize temperature to zero at all points
tt[int(N/2)] = 1./h             # Initial cond. is delta function in center
## The boundary conditions are tt[0] = tt[N-1] = 0

#* Set up loop and plot variables.
xplot = np.arange(N)*h - L/2.    # Record the x scale for plots
iplot = 0                        # Counter used to count plots
nstep = 300                      # Maximum number of iterations
nplots = 50                      # Number of snapshots (plots) to take
plot_step = nstep/nplots         # Number of time steps between plots


#* Loop over the desired number of time steps.
ttplot = np.empty((N,nplots))
tplot = np.empty(nplots)
for istep in range(nstep):  ## MAIN LOOP ##
    
    #* Compute new temperature using FTCS scheme.
    tt[1:(N-1)] = ( tt[1:(N-1)] + 
      coeff*( tt[2:N] + tt[0:(N-2)] - 2*tt[1:(N-1)] ) )
    
    #* Periodically record temperature for plotting.
    if (istep+1) % plot_step < 1 :         # Every plot_step steps
        ttplot[:,iplot] = np.copy(tt)      # record tt(i) for plotting
        tplot[iplot] = (istep+1)*tau       # Record time for plots
        iplot += 1


#* Plot temperature versus x and t as a wire-mesh plot

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection = '3d')
Tp, Xp = np.meshgrid(tplot, xplot)
ax.plot_surface(Tp, Xp, ttplot, rstride=2, cstride=2, cmap=cm.gray)
ax.set_xlabel('Time')
ax.set_ylabel('x')
ax.set_zlabel('T(x,t)')
ax.set_title('Diffusion of a delta spike')
plt.show()

#* Plot temperature versus x and t as a contour plot

levels = np.linspace(0., 10., num=21) 
ct = plt.contour(tplot, xplot, ttplot, levels) 
plt.clabel(ct, fmt='%1.2f') 
plt.xlabel('Time')
plt.ylabel('x')
plt.title('Temperature contour plot')
plt.show()
