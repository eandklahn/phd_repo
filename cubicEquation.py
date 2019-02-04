from projects import DyIDyII
import numpy as np

def getCubicCoefficients(X):
    """
    Calculates the coefficients of the characteristic polynomium for a second-order tensor.
    
    Arguments:
        X, a (3,3) Numpy array or matrix
    
    Returns:
        (a,b,c,d), where 0 = a*x**3 + b*x**2 + c*x + d
    """
    
    X1 = X[0,0]
    X2 = X[1,1]
    X3 = X[2,2]
    X4 = X[0,1]
    X5 = X[0,2]
    X6 = X[1,2]
    
    a = -1
    b = X1 + X2 + X3
    c = X4**2 + X5**2 + X6**2 - (X1*X2 + X1*X3 + X2*X3)
    d = X1*X2*X3 + 2*X4*X5*X6 - (X1*X6**2 + X2*X5**2 + X3*X4**2)
    
    return (a,b,c,d)
    
def solveCubicEquation(coeffs):
    """
    Calculates the solutions to the characteristic polynomium for a second-order tensor.
    
    Arguments:
        coeffs, (a,b,c,d) where 0 = a*x**3 + b*x**2 + c*x + d
     
    Returns:
        xes, a tuple of length 3 with solutions to the cubic equation.
    """

    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]
    d = coeffs[3]
    
    D0 = b**2 - 3*a*c
    D1 = 2*b**3 - 9*a*b*c + 27*a**2*d
    D = (D1**2 - 4*D0**3)/(-1*27*a**2)
    assert D>0 #To check that all eigenvalues will in fact be real
    C = np.power(((np.sqrt((D1**2 - 4*D0**3) + 0j) + D1)/2),1/3)
    zeta = -1/2 + 1j/2*3**(1/2)
    
    xes = []
    for n in range(3):
        x = -1/(3*a)*(b+zeta**n*C+D0/(zeta**n*C))
        xes.append(x)
    
    return tuple(xes)
    
X = DyIDyII.DyI_X_ijk
t = getCubicCoefficients(X)
p = solveCubicEquation(t)


print(p)