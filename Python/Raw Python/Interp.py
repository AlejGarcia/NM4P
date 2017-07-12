# interp - Program to interpolate data using Lagrange 
# polynomial to fit quadratic to three data points

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

def intrpf(xi,x,y):
    """Function to interpolate between data points
       using Lagrange polynomial (quadratic)
       Inputs
        x    Vector of x coordinates of data points (3 values)
        y    Vector of y coordinates of data points (3 values)
        xi   The x value where interpolation is computed
      Output
        yi   The interpolation polynomial evaluated at xi
    """

    #* Calculate yi = p(xi) using Lagrange polynomial
    yi = ( (xi-x[1])*(xi-x[2])/((x[0]-x[1])*(x[0]-x[2])) * y[0]
    + (xi-x[0])*(xi-x[2])/((x[1]-x[0])*(x[1]-x[2])) * y[1]
    + (xi-x[0])*(xi-x[1])/((x[2]-x[0])*(x[2]-x[1])) * y[2] )
    return yi

#* Initialize the data points to be fit by quadratic
x = np.empty(3)
y = np.empty(3)
print 'Enter data points as x,y pairs (e.g., [1, 2])'
for i in range(3):
    temp = np.array(input('Enter data point: '))
    x[i] = temp[0]
    y[i] = temp[1]

#* Establish the range of interpolation (from x_min to x_max)
xr = np.array(input('Enter range of x values as [x_min, x_max]: '))

#* Find yi for the desired interpolation values xi using
#  the function intrpf
nplot = 100     # Number of points for interpolation curve
xi = np.empty(nplot)
yi = np.empty(nplot)
for i in range(nplot) :
    xi[i] = xr[0] + (xr[1]-xr[0])* i/float(nplot)
    yi[i] = intrpf(xi[i], x, y)    # Use intrpf function to interpolate

#* Plot the curve given by (xi,yi) and mark original data points
plt.plot(x,y,'*',xi,yi,'-')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Three point interpolation')
plt.legend(['Data points','Interpolation  '])
plt.show()

