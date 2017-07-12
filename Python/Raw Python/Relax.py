# relax - Program to solve the Laplace equation using 
# Jacobi, Gauss-Seidel and SOR methods on a square grid

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Initialize parameters (system size, grid spacing, etc.)
method = input('Choose numerical method, 1) Jacobi; 2) Gauss-Seidel; 3) SOR')
N = input('Enter number of grid points on a side: ')
L = 1.            # System size (length)
h = L/(N-1)       # Grid spacing
x = np.arange(N)*h   # x coordinate
y = np.arange(N)*h   # y coordinate

#* Select over-relaxation factor (SOR only)
if method == 3 :
    omegaOpt = 2./(1.+np.sin(np.pi/N))    # Theoretical optimum
    print 'Theoretical optimum omega = ', omegaOpt
    omega = input('Enter desired omega: ')

#* Set initial guess as first term in separation of variables soln.
phi0 = 1.     # Potential at y=L
phi = np.empty((N,N))
for i in range(N) :
    for j in range(N) :
        phi[i,j] = phi0 * 4/(np.pi*np.sinh(np.pi)
                            ) * np.sin(np.pi*x[i]/L)*np.sinh(np.pi*y[j]/L)

#* Set boundary conditions
phi[0,:] = 0.
phi[-1,:] = 0.
phi[:,0] = 0.
phi[:,-1] = phi0*np.ones(N)    
print 'Potential at y=L equals ', phi0
print 'Potential is zero on all other boundaries'

#* Loop until desired fractional change per iteration is obtained
newphi = np.copy(phi)    # Copy of the solution (used only by Jacobi)
iterMax = N**2           # Set max to avoid excessively long runs
change = np.empty(iterMax)
changeDesired = 1.e-4    # Stop when the change is given fraction
print 'Desired fractional change = ', changeDesired
for iter in range(iterMax) :
    changeSum = 0
  
    if method == 1 :      ## Jacobi method ##
        for i in range(1,N-1) :     # Loop over interior points only
            for j in range(1,N-1) :     
                newphi[i,j] = .25*( phi[i+1,j] + phi[i-1,j] + 
                                    phi[i,j-1] + phi[i,j+1] )
                changeSum += abs( 1 - phi[i,j]/newphi[i,j] )
        phi = np.copy(newphi)   

    elif method == 2 :    ## G-S method ##
        for i in range(1,N-1) :     # Loop over interior points only
            for j in range(1,N-1) :     
                temp = .25*( phi[i+1,j] + phi[i-1,j] + 
                             phi[i,j-1] + phi[i,j+1] )
                changeSum += abs( 1 - phi[i,j]/temp )
                phi[i,j] = temp

    else :                ## SOR method ##    
        for i in range(1,N-1) :     # Loop over interior points only
            for j in range(1,N-1) :     
                temp = .25*omega*( phi[i+1,j] + phi[i-1,j] + 
                                   phi[i,j-1] + phi[i,j+1] ) + (1-omega)*phi[i,j]
                                       
                changeSum += abs( 1 - phi[i,j]/temp )
                phi[i,j] = temp

    #* Check if fractional change is small enough to halt the iteration
    change[iter] = changeSum/(N-2)**2
    if (iter+1) % 10 < 1 :
        print 'After ', iter+1, ' iterations, fractional change = ', change[iter]

    if change[iter] < changeDesired : 
        print 'Desired accuracy achieved after ', iter+1, ' iterations' 
        print 'Breaking out of main loop'
        break


#* Plot final estimate of potential as a contour plot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

levels = np.linspace(0., 1., num=11) 
ct = plt.contour(x, y, np.flipud(np.rot90(phi)), levels) 
plt.clabel(ct, fmt='%1.2f') 
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potential after %g iterations' % iter)
plt.show()

#* Plot final estimate of potential as contour and surface plots

fig = plt.figure()
ax = fig.gca(projection = '3d')
Xp, Yp = np.meshgrid(x, y)
ax.plot_surface(Xp, Yp, np.flipud(np.rot90(phi)), rstride=1, cstride=1, cmap=cm.gray)
ax.view_init(elev=30., azim=210.)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel(r'$\Phi(x,y)$')
plt.show()
