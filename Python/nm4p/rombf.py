import numpy as np

def rombf(a,b,N,func,param) :
    """Function to compute integrals by Romberg algorithm
       R = rombf(a,b,N,func,param)
       Inputs
         a,b    Lower and upper bound of the integral
         N      Romberg table is N by N
         func   Name of integrand function in a string such as
                func='errintg'.  The calling sequence is func(x,param)
         param  Set of parameters to be passed to function
       Output 
          R     Romberg table; Entry R(N,N) is best estimate of
                the value of the integral
    """
    
    #* Compute the first term R(1,1)
    h = b - a         # This is the coarsest panel size
    npanels = 1       # Current number of panels
    R = np.zeros((N+1,N+1))
    R[1,1] = h/2. * (func(a,param) + func(b,param))

    #* Loop over the desired number of rows, i = 2,...,N
    for i in range(2,N+1) :

        #* Compute the summation in the recursive trapezoidal rule
        h = h/2.          # Use panels half the previous size
        npanels *= 2      # Use twice as many panels
        sumT = 0.
        # This for loop goes k=1,3,5,...,npanels-1
        for k in range(1,npanels,2) :  
            sumT += func(a + k*h, param)
  
        #* Compute Romberg table entries R(i,1), R(i,2), ..., R(i,i)
        R[i,1] = 0.5 * R[i-1,1] + h * sumT   
        m = 1
        for j in range(2,i+1) :
            m *= 4
            R[i,j] = R[i,j-1] + ( R[i,j-1] - R[i-1,j-1] )/(m-1)

    return R