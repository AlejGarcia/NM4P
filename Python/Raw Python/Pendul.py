# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

#* Select the numerical method to use: Euler or Verlet
NumericalMethod = input('Choose a numerical method (1: Euler; 2: Verlet): ')

#* Set initial position and velocity of pendulum
theta0 = input('Enter initial angle (in degrees): ')
theta = theta0 * np.pi /180     # Convert angle to radians
omega = 0.0                     # Set the initial velocity

#* Set the physical constants and other variables
g_over_L = 1.0            # The constant g/L
time = 0.0                # Initial time
irev = 0                  # Used to count number of reversals
tau = input('Enter time step: ')

#* Take one backward step to start Verlet
accel = -g_over_L * np.sin(theta)    # Gravitational acceleration
theta_old = theta - omega*tau + 0.5*accel*tau**2    

#* Loop over desired number of steps with given time step
#    and numerical method
nstep = input('Enter number of time steps: ')
t_plot = np.empty(nstep)
th_plot = np.empty(nstep)
period = np.empty(nstep)   # Used to record period estimates
for istep in range(nstep):  

    #* Record angle and time for plotting
    t_plot[istep] = time            
    th_plot[istep] = theta * 180 / np.pi  # Convert angle to degrees
    time = time + tau
  
    #* Compute new position and velocity using 
    #    Euler or Verlet method
    accel = -g_over_L * np.sin(theta)   # Gravitational acceleration
    if NumericalMethod == 1 :
        theta_old = theta               # Save previous angle
        theta = theta + tau*omega       # Euler method
        omega = omega + tau*accel 
    else:  
        theta_new = 2*theta - theta_old + tau**2 * accel
        theta_old = theta               # Verlet method
        theta = theta_new  
  
    #* Test if the pendulum has passed through theta = 0;
    #    if yes, use time to estimate period
    if theta*theta_old < 0 :  # Test position for sign change
        print 'Turning point at time t = ',time
        if irev == 0 :          # If this is the first change,
            time_old = time     # just record the time
        else:
            period[irev-1] = 2*(time - time_old)
            time_old = time
        irev = irev + 1     # Increment the number of reversals

# Estimate period of oscillation, including error bar
nPeriod = irev-1    # Number of times the period was measured
AvePeriod = np.mean( period[0:nPeriod] )
ErrorBar = np.std(period[0:nPeriod]) / np.sqrt(nPeriod)
print 'Average period = ', AvePeriod, ' +/- ', ErrorBar

# Graph the oscillations as theta versus time
plt.plot(t_plot, th_plot, '+')
plt.xlabel('Time')
plt.ylabel(r'$\theta$ (degrees)')
plt.show()

