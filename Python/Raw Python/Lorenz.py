# lorenz - Program to compute the trajectories of the Lorenz 
# equations using the adaptive Runge-Kutta method.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

from nm4p.rk4 import rk4
from nm4p.rka import rka

#* Define the lorzrk function used by the Runge-Kutta routines
def lorzrk(s,t,param):
    """Returns right-hand side of Lorenz model ODEs
       Inputs
         s      State vector [x y z]
         t      Time (not used)
         param  Parameters [r sigma b]
       Output
         deriv  Derivatives [dx/dt dy/dt dz/dt]
    """
    
    #* For clarity, unravel input vectors
    x, y, z = s[0], s[1], s[2]
    r = param[0]
    sigma = param[1]
    b = param[2]

    #* Return the derivatives [dx/dt dy/dt dz/dt]
    deriv = np.empty(3)
    deriv[0] = sigma*(y-x)
    deriv[1] = r*x - y - x*z
    deriv[2] = x*y - b*z
    return deriv


#* Set initial state x,y,z and parameters r,sigma,b
state = np.array(input('Enter the initial position [x, y, z]: '))
r = input('Enter the parameter r: ')
sigma = 10.    # Parameter sigma
b = 8./3.      # Parameter b
param = np.array([r, sigma, b])  # Vector of parameters passed to rka
tau = 1.       # Initial guess for the timestep
err = 1.e-3    # Error tolerance

#* Loop over the desired number of steps
time = 0.
nstep = input('Enter number of steps: ')
tplot = np.empty(nstep)
tauplot = np.empty(nstep)
xplot, yplot, zplot = np.empty(nstep), np.empty(nstep), np.empty(nstep)
for istep in range(nstep):

    #* Record values for plotting
    x, y, z = state[0], state[1], state[2]
    tplot[istep] = time
    tauplot[istep] = tau       
    xplot[istep] = x;    yplot[istep] = y;    zplot[istep] = z 
    if (istep+1) % 50  < 1 :
        print 'Finished ',istep, ' steps out of ',nstep

    #* Find new state using adaptive Runge-Kutta
    [state, time, tau] = rka(state, time, tau, err, lorzrk, param)


#* Print max and min time step returned by rka
tauMax = np.max(tauplot[1:nstep])
tauMin = np.min(tauplot[1:nstep])
print 'Adaptive time step: Max = ', tauMax, ' Min = ', tauMin

#* Graph the time series x(t)
plt.plot(tplot,xplot,'-')
plt.xlabel('Time')
plt.ylabel('x(t)')
plt.title('Lorenz model time series')
plt.show()

#* Graph the x,y,z phase space trajectory
# Mark the location of the three steady states
x_ss = np.empty(3);  y_ss = np.empty(3);  z_ss = np.empty(3)
x_ss[0] = 0.
y_ss[0] = 0.
z_ss[0] = 0.
x_ss[1] = np.sqrt( b*(r-1.) )
y_ss[1] = x_ss[1]
z_ss[1] = r - 1.
x_ss[2] = -np.sqrt( b*(r-1.) )
y_ss[2] = x_ss[2]
z_ss[2] = r - 1.

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(xplot, yplot, zplot ,'-')
ax.plot(x_ss, y_ss, z_ss, '*')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Lorenz model phase space')
plt.show()

