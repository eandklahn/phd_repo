import numpy as np
import scientificConstants as SC
import matplotlib.pyplot as plt

def effHalfSpinG(X, sX, T):
    """See the document 'g-tensor from CHI in the effective spin½-model'
    sX can be determined with the sampling method
    
    Input
    X: Easy-axis magnetic susceptibility
    sX: Uncertainty on the easy-axis magnetic susceptibility
    T: absolute temperature
    
    Output
    g: easy-axis g value
    sg: uncertainty on the easy-axis g-value
    """
    
    g = np.sqrt(4*SC.kB*X*T/SC.muB)
    
    sg = np.sqrt(SC.kB*T/(SC.muB*X))*sX
    
    return (g,sg)

def coth(x):
    """Creates the hyperbolic cotangent function from Numpys implementation of sinh and cosh"""
    r = np.cosh(x)/np.sinh(x)
    return r

def landeG(S,L,J):
    """Function that returns the Landé g-factor for a tuple of (S,L,J)"""
    g = 3/2 + (S*(S+1)-L*(L+1))/(2*J*(J+1))
    return g

def brillouin(J, y):
    """Returns the Brillouin function of J as a function of y"""
    
    assert type(y) == np.ndarray
    assert type(J) == float
    assert y.ndim == 1
    
    frac = (2*J+1)/(2*J)
    BRILL = frac*coth(frac*y)-(1/(2*J))*coth(y/(2*J))
    
    return BRILL
    
def brillouinPlot(J, L, S, maxB, maxM, T, minB=1e-6, increments=100):
    """Function that returns a tuple containing B and y (for x-axis) and M (for y-axis) for a specific maximum
    magnetization, using the equations found on p. 28 of Blundell.
    J: total angular momentum of state
    L: orbital angular momentum of state
    S: spin angular momentum of state
    maxB: the maximum magnetic field to plot after
    maxM: saturation magnetization
    T: temperature (in Kelvin)"""
    
    #An array of magnetic field values (unit: Tesla)
    B = np.linspace(minB, maxB, increments)
    
    #The array to use in the Brillouin function
    y = landeG(J, S, L)*SC.muB*J/(SC.kB*T)*B
    
    #Magnetization array calculated with the saturation magnetization and the Brillouin function
    M = maxM*brillouin(J,y)
    
    return (y,B,M)





