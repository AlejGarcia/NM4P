import numpy as np

def pollsf(x, y, sigma, M):
    """Function to fit a polynomial to data
       Inputs 
        x       Independent variable
        y       Dependent variable
        sigma   Estimate error in y
        M       Number of parameters used to fit data
       Outputs
        a_fit   Fit parameters; a(1) is intercept, a(2) is slope
        sig_a   Estimated error in the parameters a()
        yy      Curve fit to the data
        chisqr  Chi squared statistic
    """
    
    #* Form the vector b and design matrix A   
    N = len(x)
    b = np.empty(N)
    A = np.empty((N,M))
    for i in range(N):
        b[i] = y[i]/sigma[i]
        for j in range(M):
            A[i,j] = x[i]**j / sigma[i] 

    #* Compute the correlation matrix C 
    C = np.linalg.inv( np.dot( np.transpose(A), A) )

    #* Compute the least squares polynomial coefficients a_fit
    a_fit = np.dot(C, np.dot( np.transpose(A), np.transpose(b)) )

    #* Compute the estimated error bars for the coefficients
    sig_a = np.empty(M)
    for j in range(M):
        sig_a[j] = np.sqrt(C[j,j])

    #* Evaluate curve fit at each data point and compute Chi^2
    yy = np.zeros(N)
    chisqr = 0.
    for i in range(N):
        for j in range(M):
            yy[i] += a_fit[j]*x[i]**j   # yy is the curve fit
        chisqr += ((y[i]-yy[i]) / sigma[i])**2
        
    return [a_fit, sig_a, yy, chisqr]