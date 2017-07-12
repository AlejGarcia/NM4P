import numpy as np

def tri_ge(a,b) :
    # Function to solve b = a*x by Gaussian elimination where
    # the matrix a is a packed tridiagonal matrix
    # Inputs
    #   a    Packed tridiagonal matrix, N by N unpacked
    #   b    Column vector of length N
    # Output 
    #   x    Solution of b = a*x; Column vector of length N

    #* Check that dimensions of a and b are compatible
    N_a = np.shape(a)
    N = len(b)
    if N_a[0] != N or N_a[1] != 3 :
        print 'Problem in tri_GE, inputs are incompatible'
        return None

    #* Unpack diagonals of triangular matrix into vectors
    alpha = np.copy(a[1:N,0])
    beta = np.copy(a[:,1])
    gamma = np.copy(a[0:(N-1),2])
    bb = np.copy(b)  # Copy is modified in forward elimination

    #* Perform forward elimination
    for i in range(1,N) :
        coeff = alpha[i-1]/beta[i-1]
        beta[i] -= coeff*gamma[i-1]
        bb[i] -= coeff*bb[i-1]

    #* Perform back substitution
    x = np.empty(N,dtype=complex)
    x[-1] = bb[-1]/beta[-1]
    for i in reversed(range(N-1)) :
        x[i] = (bb[i] - gamma[i] * x[i+1])/beta[i]

    return x