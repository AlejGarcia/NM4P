# neutrn - Program to solve the neutron diffusion equation 
# using the Forward Time Centered Space (FTCS) scheme.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Initialize parameters (time step, grid points, etc.).
tau = input('Enter time step: ')
N = input('Enter the number of grid points: ')
L = input('Enter system length: ')
# The system extends from x=-L/2 to x=L/2
h = L/float(N-1)   # Grid size
D = 1.    # Diffusion coefficient
C = 1.    # Generation rate
coeff = D*tau/h**2
coeff2 = C*tau 
if coeff < 0.5 :
    print 'Solution is expected to be stable'
else:
    print 'WARNING: Solution is expected to be unstable'

#* Set initial and boundary conditions.
nn = np.zeros(N)           # Initialize density to zero at all points
nn_new = np.zeros(N)       # Initialize temporary array used by FTCS
nn[int(N/2.)] = 1/h        # Initial cond. is delta function in center
## The boundary conditions are nn[0] = nn[N-1] = 0

#* Set up loop and plot variables.
xplot = np.arange(N)*h - L/2.    # Record the x scale for plots
iplot = 0                        # Counter used to count plots
nstep = input('Enter number of time steps: ')
nplots = 50                # Number of snapshots (plots) to take
plot_step = nstep/nplots   # Number of time steps between plots

#* Loop over the desired number of time steps.
nnplot = np.empty((N,nplots))
tplot = np.empty(nplots)
nAve = np.empty(nplots)
for istep in range(nstep):     ## MAIN LOOP ##

    #* Compute the new density using FTCS scheme.
    nn[1:(N-1)] = ( nn[1:(N-1)] + 
        coeff*( nn[2:N] + nn[0:(N-2)] - 2*nn[1:(N-1)] ) +
        coeff2*nn[1:(N-1)] )
  
    #* Periodically record the density for plotting.
    if (istep+1) % plot_step < 1:      # Every plot_step steps
        nnplot[:,iplot] = np.copy(nn)  # record nn[i] for plotting
        tplot[iplot] = (istep+1)*tau   # Record time for plots
        nAve[iplot] = np.mean(nn)      # Record average density 
        iplot += 1 
        print 'Finished ', istep, ' of ', nstep, ' steps'


#* Plot density versus x and t as a 3D-surface plot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection = '3d')
Tp, Xp = np.meshgrid(tplot, xplot)
ax.plot_surface(Tp, Xp, nnplot, rstride=2, cstride=2, cmap=cm.gray)
ax.set_xlabel('Time')
ax.set_ylabel('x')
ax.set_zlabel('n(x,t)');
ax.set_title('Neutron diffusion');
plt.show()

#* Plot average neutron density versus time
plt.plot(tplot,nAve,'*')
plt.xlabel('Time')
plt.ylabel('Average density')
plt.title(r'$L$ = %g,  ($L_c = \pi$)' % L)
plt.show()

