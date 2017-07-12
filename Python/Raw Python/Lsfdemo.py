# lsfdemo - Program for demonstrating least squares fit routines

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

from nm4p.linreg import linreg
from nm4p.pollsf import pollsf

#* Initialize data to be fit. Data is quadratic plus random number.
print 'Curve fit data is created using the quadratic'
print '  y(x) = c(0) + c(1)*x + c(2)*x**2'
c = np.array(input('Enter the coefficients as [c(0), c(1), c(2)]: '))
N = 50;                 # Number of data points
x = np.arange(1,N+1)    # x = [1, 2, ..., N]
y = np.empty(N)
alpha = input('Enter estimated error bar: ')
sigma = alpha * np.ones(N)  # Constant error bar
np.random.seed(0)           # Initialize random number generator
for i in range(N):
    r = alpha * np.random.normal()    # Gaussian distributed random vector
    y[i] = c[0] + c[1]*x[i] + c[2]*x[i]**2 + r       

#* Fit the data to a straight line or a more general polynomial
M = input('Enter number of fit parameters (=2 for line): ')
if M == 2 :  
    #* Linear regression (Straight line) fit
    [a_fit, sig_a, yy, chisqr] = linreg(x, y, sigma)
else: 
    #* Polynomial fit
    [a_fit, sig_a, yy, chisqr] = pollsf(x, y, sigma, M)


#* Print out the fit parameters, including their error bars.
print 'Fit parameters:'
for i in range(M):
    print ' a[', i, '] = ', a_fit[i], ' +/- ', sig_a[i]

#* Graph the data, with error bars, and fitting function.
plt.errorbar(x,y,sigma,None,'o')   # Graph data with error bars
plt.plot(x,yy,'-')            # Plot the fit on same graph as data
plt.xlabel(r'$x_i$')  
plt.ylabel(r'$y_i$ and $Y(x)$') 
plt.title('chi^2 = %d,    N-M = %d' % (chisqr, N-M) )
plt.show()

