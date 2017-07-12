# newtn - Program to solve a system of nonlinear equations 
# using Newton's method.  Equations defined by function fnewt.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

# Define the function f(x,a) used for root finding
def fnewt(x,a):
    """Function used by the N-variable Newton's method
       Inputs
         x     State vector [x y z]
         a     Parameters [r sigma b]
       Outputs
         f     Lorenz model r.h.s. [dx/dt dy/dt dz/dt]
         D     Jacobian matrix, D(i,j) = df(j)/dx(i)
    """
    
    # Evaluate f(i)
    f = np.empty(3)
    f[0] = a[1] * (x[1]-x[0])
    f[1] = a[0]*x[0] -x[1] -x[0]*x[2]
    f[2] = x[0]*x[1] -a[2]*x[2]

    # Evaluate D(i,j)
    D = np.empty((3,3))
    D[0,0] = -a[1]        # df(0)/dx(0)
    D[0,1] = a[0]-x[2]    # df(1)/dx(0)
    D[0,2] = x[1]         # df(2)/dx(0)
    D[1,0] = a[1]         # df(0)/dx(1)
    D[1,1] = -1.          # df(1)/dx(1)
    D[1,2] = x[0]         # df(2)/dx(1)
    D[2,0] = 0.           # df(0)/dx(2)
    D[2,1] = -x[0]        # df(1)/dx(2)
    D[2,2] = -a[2]        # df(2)/dx(2)

    return [f, D]


#* Set initial guess and parameters
x0 = np.array(input('Enter the initial guess (row vector): '))
x = np.copy(x0)     # Copy initial guess
a = np.array(input('Enter the parameter a: '))

#* Loop over desired number of steps 
nStep = 10    # Number of iterations before stopping
xp = np.empty((len(x), nStep))
xp[:,0] = np.copy(x[:])     # Record initial guess for plotting
for iStep in range(nStep):

    #* Evaluate function f and its Jacobian matrix D
    [f, D] = fnewt(x,a)       # fnewt returns value of f and D
    
    #* Find dx by Gaussian elimination; transpose D for column vectors
    dx = np.linalg.solve( np.transpose(D), f)    
    
    #* Update the estimate for the root  
    x = x - dx                    # Newton iteration for new x
    xp[:,iStep] = np.copy(x[:])   # Save current estimate for plotting


#* Print the final estimate for the root
print 'After', nStep, ' iterations the root is'
print x

# %* Plot the iterations from initial guess to final estimate
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(xp[0,:],xp[1,:],xp[2,:],'*-')
ax.plot([x[0]],[x[1]],[x[2]],'ro')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Steady state of the Lorenz model')
plt.show()
