# sprfft - Program to compute the power spectrum of a  
# coupled mass-spring system.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

from nm4p.rk4 import rk4

# Define the sprrk function used by the Runge-Kutta routines
def sprrk(s,t,param):
    """Returns right-hand side of 3 mass-spring system
       equations of motion
       Inputs
         s       State vector [x(1) x(2) ... v(3)]
         t       Time (not used)
         param   (Spring constant)/(Block mass)
       Output
         deriv   [dx(1)/dt dx(2)/dt ... dv(3)/dt]
    """
    
    deriv = np.empty(6)
    deriv[0] = s[3]
    deriv[1] = s[4]
    deriv[2] = s[5]
    param2 = -2.*param
    deriv[3] = param2*s[0] + param*s[1]
    deriv[4] = param2*s[1] + param*(s[0]+s[2])
    deriv[5] = param2*s[2] + param*s[1]
    return deriv


#* Set parameters for the system (initial positions, etc.).
x = np.array(input('Enter initial displacement [x0, x1, x2]: '))  
v = np.array([0., 0., 0.])       # Masses are initially at rest
# Positions and velocities; used by rk4
state = np.array([x[0], x[1], x[2], v[0], v[1], v[2]])      
tau = input('Enter timestep: ')  
k_over_m = 1.               # Ratio of spring const. over mass

#* Loop over the desired number of time steps.
time = 0.          # Set initial time
nstep = 256        # Number of steps in the main loop
nprint = nstep/8   # Number of steps between printing progress
tplot = np.empty(nstep)
xplot = np.empty((nstep,3))
for istep in range(nstep):  ### MAIN LOOP ###

    #* Use Runge-Kutta to find new displacements of the masses.
    state = rk4(state,time,tau,sprrk,k_over_m)  
    time = time + tau
  
    #* Record the positions for graphing and to compute spectra.
    xplot[istep,:] = np.copy(state[0:3])   # Record positions
    tplot[istep] = time
    if istep % nprint < 1 :
        print 'Finished ', istep, ' out of ', nstep, ' steps'


#* Graph the displacements of the three masses.
plt.plot(tplot,xplot[:,0],'-',tplot,xplot[:,1],'-.',tplot,xplot[:,2],'--')
plt.legend(['Mass #1  ','Mass #2  ','Mass #3  '])
plt.title('Displacement of masses (relative to rest positions)')
plt.xlabel('Time') 
plt.ylabel('Displacement')
plt.show()

#* Calculate the power spectrum of the time series for mass #1
f = np.arange(nstep)/(tau*nstep)   # Frequency
x1 = xplot[:,0]                # Displacement of mass 1

x1fft = np.fft.fft(x1)       # Fourier transform of displacement

spect = np.empty(len(x1fft))        # Power spectrum of displacement
for i in range(len(x1fft)):
    spect[i] = abs(x1fft[i])**2


#* Apply the Hanning window to the time series and calculate
#  the resulting power spectrum
x1w = np.empty(len(x1))
for i in range(len(x1)):
    window = 0.5 * (1. - np.cos(2*np.pi*i/nstep)) # Hanning window
    x1w[i] = x1[i] * window          # Windowed time series
    
x1wfft = np.fft.fft(x1w)            # Fourier transf. (windowed data)

spectw = np.empty(len(x1wfft))      # Power spectrum (windowed data)
for i in range(len(x1wfft)):
    spectw[i] = abs(x1wfft[i])**2


#* Graph the power spectra for original and windowed data
plt.semilogy(f[0:(nstep/2)],spect[0:(nstep/2)],'-', 
             f[0:(nstep/2)],spectw[0:(nstep/2)],'--')
plt.title('Power spectrum (dashed is windowed data)')
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.show()
