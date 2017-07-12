import numpy as np

class sortList:
    """Class used for sorting particles into cells."""
    def __init__(self, ncell_in, npart_in):
        self.ncell = ncell_in
        self.npart = npart_in
        self.cell_n = np.zeros(ncell_in, dtype=int)
        self.index = np.empty(ncell_in, dtype=int)
        self.Xref = np.empty(npart_in, dtype=int)
        
class sampList:
    """Class used for sampling density, velocity, and temperature."""
    def __init__(self, ncell_in):
        self.ncell = ncell_in
        self.nsamp = 0
        self.ave_n = np.zeros(ncell_in)
        self.ave_u = np.zeros((ncell_in,3))
        self.ave_T = np.zeros(ncell_in)

#################################################################        
def sorter(x,L,sD) :
    """sorter - Function to sort particles into cells
       Inputs
         x       Positions of particles
         L       System size
         sD      Object containing sorting lists
    """
    
    #* Find the cell address for each particle
    npart = sD.npart
    ncell = sD.ncell
    jx = np.empty(npart,dtype=int)
    for ipart in range(npart) :
        jx[ipart] = int( x[ipart]*ncell/L )
        jx[ipart] = min( jx[ipart], (ncell-1) )

    #* Count the number of particles in each cell
    sD.cell_n = np.zeros(ncell)
    for ipart in range(npart) :
        sD.cell_n[ jx[ipart] ] += 1

    #* Build index list as cumulative sum of the 
    #  number of particles in each cell
    m = 0
    for jcell in range(ncell) :
        sD.index[jcell] = m
        m += sD.cell_n[jcell]

    #* Build cross-reference list
    temp = np.zeros(ncell, dtype=int)      # Temporary array
    for ipart in range(npart) :
        jcell = jx[ipart]       # Cell address of ipart
        k = sD.index[jcell] + temp[jcell]
        sD.Xref[k] = ipart
        temp[jcell] += 1

#################################################################
def colider(v,crmax,tau,selxtra,coeff,sD) :
    """colider - Function to process collisions in cells
       Inputs
         v         Velocities of the particles
         crmax     Estimated maximum relative speed in a cell
         tau       Time step
         selxtra   Extra selections carried over from last timestep
         coeff     Coefficient in computing number of selected pairs
         sD        Object containing sorting lists 
       Outputs
         col       Total number of collisions processed
    """
    
    ncell = sD.ncell 
    col = 0              # Count number of collisions
    vrel = np.empty(3)   # Relative velocity for collision pair
    
    #* Loop over cells, processing collisions in each cell
    for jcell in range(ncell) :
 
        #* Skip cells with only one particle
        number = sD.cell_n[jcell]
        if number > 1 :  
            
            #* Determine number of candidate collision pairs 
            #  to be selected in this cell
            select = coeff*number*(number-1)*crmax[jcell] + selxtra[jcell]
            nsel = int(select)            # Number of pairs to be selected
            selxtra[jcell] = select-nsel  # Carry over any left-over fraction
            crm = crmax[jcell]            # Current maximum relative speed
  
            #* Loop over total number of candidate collision pairs
            for isel in range(nsel) :
    
                #* Pick two particles at random out of this cell
                k = int( np.floor( np.random.uniform(0,number) ) )
                kk = int(np.ceil( k + np.random.uniform(0,number-1) ) % number )
                ip1 = sD.Xref[ k + sD.index[jcell] ]   # First particle
                ip2 = sD.Xref[ kk + sD.index[jcell] ]  # Second particle

                #* Calculate pair's relative speed
                cr = np.linalg.norm( v[ip1,:] - v[ip2,:] )   # Relative speed 
                if cr > crm :         # If relative speed larger than crm,
                    crm = cr          # then reset crm to larger value

                #* Accept or reject candidate pair according to relative speed
                if cr/crmax[jcell] > np.random.random()  :
                    #* If pair accepted, select post-collision velocities
                    col += 1                            # Collision counter
                    vcm = 0.5*( v[ip1,:] + v[ip2,:] )   # Center of mass velocity
                    cos_th = 1. - 2.*np.random.random()    # Cosine and sine of 
                    sin_th = np.sqrt(1. - cos_th**2)    # collision angle theta
                    phi = 2*np.pi*np.random.random()       # Collision angle phi
                    vrel[0] = cr*cos_th                 # Compute post-collision 
                    vrel[1] = cr*sin_th*np.cos(phi)     # relative velocity
                    vrel[2] = cr*sin_th*np.sin(phi)
                    v[ip1,:] = vcm + 0.5*vrel           # Update post-collision
                    v[ip2,:] = vcm - 0.5*vrel           # velocities

            crmax[jcell] = crm      # Update max relative speed 
    
    return col

#################################################################
def mover( x, v, npart, L, mpv, vwall, tau) :
    """mover - Function to move particles by free flight
               Also handles collisions with walls
       Inputs
         x        Positions of the particles
         v        Velocities of the particles
         npart    Number of particles in the system
         L        System length
         mpv      Most probable velocity off the wall
         vwall    Wall velocities
         tau      Time step
       Outputs
         strikes  Number of particles striking each wall
         delv     Change of y-velocity at each wall     
    """
    
    #* Move all particles pretending walls are absent
    x_old = np.copy(x)     # Remember original position
    x[:] = x_old[:] + v[:,0]*tau   

    #* Loop over all particles
    strikes = np.array([0, 0])
    delv = np.array([0., 0.])  
    xwall = np.array([0., L])
    vw = np.array([-vwall, vwall])
    direction = [1, -1]   # Direction of particle leaving wall
    stdev = mpv/np.sqrt(2.)
    for i in range(npart) :

        #* Test if particle strikes either wall
        if x[i] <= 0. :
            flag = 0   # Particle strikes left wall
        elif x[i] >= L :
            flag = 1   # Particle strikes right wall
        else :
            flag = -1   # Particle strikes neither wall

        #* If particle strikes a wall, reset its position
        #  and velocity. Record velocity change.
        if flag > -1 :
            strikes[flag] += 1
            vyInitial = v[i,1]
            #* Reset velocity components as biased Maxwellian,
            #  Exponential dist. in x; Gaussian in y and z
            v[i,0] = direction[flag] * np.sqrt(
                -np.log( 1. - np.random.random() ) ) * mpv
            v[i,1] = stdev*np.random.normal() + vw[flag]  # Add wall velocity
            v[i,2] = stdev*np.random.normal()
            # Time of flight after leaving wall
            dtr = tau * (x[i] - xwall[flag])/(x[i] - x_old[i])   
            #* Reset position after leaving wall
            x[i] = xwall[flag] + v[i,0]*dtr
            #* Record velocity change for force measurement
            delv[flag] += v[i,1] - vyInitial
            
    return [ strikes, delv ]

#################################################################
def sampler(x,v,npart,L,sampD) :
    """ sampler - Function to sample density, velocity and temperature
      Inputs
         x       Particle positions
         v       Particle velocities
         npart   Number of particles
         L       System size
         sampD   Structure with sampling data
    """
    
    #* Compute cell location for each particle
    ncell = sampD.ncell
    jx = np.empty(npart, dtype=int)
    for i in range(npart) :
        jx[i] = int(ncell*x[i]/L)

    #* Initialize running sums of number, velocity and v^2
    sum_n = np.zeros(ncell)
    sum_v = np.zeros((ncell,3))
    sum_v2 = np.zeros(ncell)

    #* For each particle, accumulate running sums for its cell
    for ipart in range(npart) :
        jcell = jx[ipart]    # Particle ipart is in cell jcell
        sum_n[jcell] += 1.
        sum_v[jcell,:] += v[ipart,:]
        sum_v2[jcell] += v[ipart,0]**2 + v[ipart,1]**2 + v[ipart,2]**2

    #* Use current sums to update sample number, velocity 
    #  and temperature
    for i in range(3) :
        sum_v[:,i] /= sum_n[:]
    sum_v2[:] /= sum_n[:]
    
    sampD.ave_n[:] += sum_n[:]
    for i in range(3) :    
        sampD.ave_u[:,i] += sum_v[:,i]
    sampD.ave_T[:] += sum_v2[:] - (
                   sum_v[:,0]**2 + sum_v[:,1]**2 + sum_v[:,2]**2 )
    sampD.nsamp += 1
 
