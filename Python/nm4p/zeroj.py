import numpy as np
from bess import bess

def zeroj( m_order, n_zero) :
    """Zeros of the Bessel function J(x)
       Inputs
        m_order = Order of the Bessel function
        n_zero  = Index of the zero (first, second, etc.)
       Output
        z = The "n_zero th" zero of the Bessel function
    """
    
    #* Use asymtotic formula for initial guess
    beta = (n_zero + 0.5*m_order - 0.25)*np.pi
    mu = 4*m_order**2
    z = beta - (mu-1.)/(8.*beta) - 4.*(mu-1.)*(7.*mu-31.)/(3.*(8.*beta)**3)

    #* Use Newton's method to locate the root
    jj = np.empty(m_order+2)
    for i in range(5) :
        jj = bess(m_order+1,z)     
        # Use the recursion relation to evaluate derivative
        deriv = -jj[m_order+1] + m_order/z * jj[m_order]
        z = z - jj[m_order]/deriv   # Newton's root finding  
    
    return z