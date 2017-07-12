# dsmcne - Program to simulate a dilute gas using DSMC algorithm
# This version simulates planar Couette flow

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

from nm4p.dsmc import sampList
from nm4p.dsmc import sortList
from nm4p.dsmc import sorter
from nm4p.dsmc import colider
from nm4p.dsmc import mover
from nm4p.dsmc import sampler


#* Initialize constants  (particle mass, diameter, etc.)
boltz = 1.3806e-23     # Boltzmann's constant (J/K)
mass = 6.63e-26        # Mass of argon atom (kg)
diam = 3.66e-10        # Effective diameter of argon atom (m)
T = 273.               # Temperature (K)
density = 2.685e25     # Number density of argon at STP (m^-3)
L = 1.e-6              # System size is one micron
Volume = L**3          # Volume of the system (m^3)
npart = input('Enter number of simulation particles: ')
eff_num = density*Volume/npart
print 'Each simulation particle represents ', eff_num, ' atoms'
mfp = Volume/(np.sqrt(2.)*np.pi*diam**2*npart*eff_num)
print 'System width is ', L/mfp, ' mean free paths'
mpv = np.sqrt(2*boltz*T/mass)   # Most probable initial velocity 
vwall_m = input('Enter wall velocity as Mach number: ')
vwall = vwall_m * np.sqrt(5./3. * boltz*T/mass)
print 'Wall velocities are ', -vwall, ' and ', vwall, ' m/s'

#* Assign random positions and velocities to particles
np.random.seed(0)          # Initialize random number generator
x = np.empty(npart)
for i in range(npart) :
    x[i] = np.random.uniform(0.,L)       # Assign random positions
v = np.zeros((npart,3))     
for i in range(npart) :
    for j in range(3) :
        # Assign thermal velocities using Gaussian random numbers
        v[i,j] = np.sqrt(boltz*T/mass) * np.random.normal()
    v[i,1] += 2. * vwall * x[i]/L - vwall   # Add velocity gradient

#* Initialize variables used for evaluating collisions
ncell = 20                       # Number of cells
tau = 0.2*(L/ncell)/mpv          # Set timestep tau
vrmax = 3*mpv*np.ones(ncell)     # Estimated max rel. speed in a cell
selxtra = np.zeros(ncell)        # Used by collision routine "colider"
coeff = 0.5*eff_num*np.pi*diam**2*tau/(Volume/ncell)

#* Declare sortList object for lists used in sorting
sortData = sortList(ncell, npart)

#* Initialize object and variables used in statistical sampling
sampData = sampList(ncell)
tsamp = 0.                # Total sampling time
dvtot = np.zeros(2)       # Total momentum change at a wall
dverr = np.zeros(2)       # Used to find error in dvtot

#* Loop for the desired number of time steps
colSum = 0
strikeSum = np.array([0, 0])
nstep = input('Enter total number of time steps: ')
for istep in range(nstep) :

    #* Move all the particles 
    [ strikes, delv ] = mover(x,v,npart,L,mpv,vwall,tau)
    strikeSum += strikes

    #* Sort the particles into cells
    sorter(x,L,sortData)
  
    #* Evaluate collisions among the particles
    col = colider(v,vrmax,tau,selxtra,coeff,sortData)
    colSum += col 
  
    #* After initial transient, accumulate statistical samples
    if istep > nstep/10 : 
        sampler(x,v,npart,L,sampData)
        dvtot += delv
        dverr += delv**2
        tsamp += tau

    #* Periodically display the current progress
    if (istep+1) % 100 < 1 :
        print 'Finished ', istep, ' of ', nstep, ' steps, Collisions = ',colSum
        print 'Total wall strikes: ', strikeSum[0], ' (left)  ', strikeSum[1], ' (right)'


#* Normalize the accumulated statistics
nsamp = sampData.nsamp 
ave_n = (eff_num/(Volume/ncell))*sampData.ave_n/nsamp
ave_u = np.empty((ncell,3))
for i in range(3) :
    ave_u[:,i] = sampData.ave_u[:,i]/nsamp
ave_T = mass/(3*boltz) * (sampData.ave_T/nsamp)
dverr = dverr/(nsamp-1) - (dvtot/nsamp)**2
dverr = np.sqrt(dverr*nsamp)

#* Compute viscosity from drag force on the walls
force = (eff_num*mass*dvtot)/(tsamp*L**2)
ferr = (eff_num*mass*dverr)/(tsamp *L**2)
print 'Force per unit area is'
print 'Left wall:   ', force[0], ' +/- ', ferr[0]  
print 'Right wall:  ', force[1], ' +/- ', ferr[1]  
vgrad = 2*vwall/L;  # Velocity gradient
visc = 1./2.*(-force[0]+force[1])/vgrad   # Average viscosity
viscerr = 1./2.*(ferr[0]+ferr[1])/vgrad   # Error
print 'Viscosity = ', visc, ' +/- ', viscerr, ' N s/m^2'
eta = 5.*np.pi/32.*mass*density*(2./np.sqrt(np.pi)*mpv)*mfp
print 'Theoretical value of viscoisty is ', eta, ' N s/m^2'

#* Plot average density, velocity and temperature
xcell = (np.arange(ncell)+0.5)/ncell * L
plt.plot(xcell,ave_n)           
plt.xlabel('position')  
plt.ylabel('Number density') 
plt.show()
plt.plot(xcell,ave_u[:,0],xcell,ave_u[:,1],xcell,ave_u[:,2])           
plt.xlabel('position')  
plt.ylabel('Velocities')
plt.legend(['x-component','y-component','z-component'])
plt.show()
plt.plot(xcell,ave_T)           
plt.xlabel('position')  
plt.ylabel('Temperature') 
plt.show()

