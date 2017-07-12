import numpy as np

def legndr(n,x) :
    """Legendre polynomials function
       Inputs  
         n = Highest order polynomial returned
         x = Value at which polynomial is evaluated
       Output
         p = Vector containing P(x) for order 0,1,...,n
    """
    
    #* Perform upward recursion
    p = np.empty(n+1)
    p[0] = 1.      # P(x) for n=0
    if n == 0 :
        return p
    p[1] = x       # P(x) for n=1
    if n == 1 :
        return p
    
    # Use upward recursion to obtain other n's
    for i in range(1,n) :
        p[i+1] = ((2.*i+1.)*x*p[i] - i*p[i-1])/(i+1.)

    return p