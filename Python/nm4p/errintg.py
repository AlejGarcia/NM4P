import numpy as np

def errintg(x,param) :
    """Error function integrand
       Inputs
         x       Value where integrand is evaluated
         param   Parameter list (not used)
       Output
         f       Integrand of the error function
    """
    
    f = np.exp(-x**2)
    return f