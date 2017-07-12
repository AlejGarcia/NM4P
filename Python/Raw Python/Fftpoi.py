# fftpoi - Program to solve the Poisson equation using 
# MFT method (periodic boundary conditions)

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt


#* Initialize parameters (system size, grid spacing, etc.)
eps0 = 8.8542e-12    # Permittivity (C^2/(N m^2))
N = 50    # Number of grid points on a side (square grid)
L = 1.    # System size
h = L/N   # Grid spacing for periodic boundary conditions
x = (np.arange(N) + 1./2)*h    # Coordinates  of grid points
y = np.copy(x)                 # Square grid
print 'System is a square of length ', L

#* Set up charge density rho(i,j) 
rho = np.zeros((N,N));  # Initialize charge density to zero
M = input('Enter number of line charges: ')
for i in range(M) :
    print '  For charge #', i
    r = input('Enter position [x, y]: ') 
    ii=int(r[0]/h + 0.5)    # Place charge at nearest
    jj=int(r[1]/h + 0.5)    # grid point
    q = input('Enter charge density: ') 
    rho[ii,jj] += q/h**2

#* Compute matrix P
cx = np.cos( (2*np.pi/N) * np.arange(N) )
cy = np.copy(cx)
numerator = -h**2/(2.*eps0)
tinyNumber = 1e-20;  # Avoids division by zero
P = np.empty((N,N))
for i in range(N) :
    for j in range(N) :
        P[i,j] = numerator/(cx[i]+cy[j]-2.+tinyNumber)

#* Compute potential using MFT method
rhoT = np.fft.fft2(rho)    # Transform rho into wavenumber domain
phiT = rhoT * P            # Computing phi in the wavenumber domain
phi = np.fft.ifft2(phiT);  # Inv. transf. phi into the coord. domain
phi = np.real(phi);        # Clean up imaginary part due to round-off

#* Compute electric field as E = - grad phi
[Ex, Ey] = np.gradient(np.flipud(np.rot90(phi))) 
for i in range(N) :
    for j in range(N) :
        magnitude = np.sqrt(Ex[i,j]**2 + Ey[i,j]**2)         
        Ex[i,j] /= -magnitude     # Normalize components so
        Ey[i,j] /= -magnitude     # vectors have equal length


#* Plot potential
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection = '3d')
Xp, Yp = np.meshgrid(x, y)
ax.contour(Xp,Yp,np.flipud(np.rot90(phi)),35)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel(r'$\Phi(x,y)$')
plt.show()

#* Plot electric field
plt.quiver(Xp,Yp,Ey,Ex)        # Plot E field with vectors
plt.title('E field (Direction)') 
plt.xlabel('x')
plt.ylabel('y')
plt.axis('square')  
plt.axis([0., L, 0., L])
plt.show()

