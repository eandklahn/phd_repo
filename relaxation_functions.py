import numpy as np

def fuki_density(s, a):
    """
    Implements the Fuoss-Kirkwood distribution for
    magnetic relaxation times.
    See also D. Reta & N. Chilton. DOI: 10.1039/C9CP04301B
    
    Arguments
    ----------
    s : float
    The value at which the distribution should be probed
    a : float
    alpha (distribution of relaxation times)
    
    Returns
    ----------
    p : float
    The value of the Fuoss-Kirkwood density at the specified t
    """
    
    p = np.sin(a*np.pi)/(2*np.pi*(np.cosh((1-a)*s)-np.cos(a*np.pi)))
    
    return p
    
def A(I, a):
    """
    Implements the integration limit A as a function of the
    value of the integral I.
    
    NOTE: This function is not working correctly yet
    
    Arguments
    ----------
    I : float
    The fraction of the full area that the integral has to encompass (1sigma, 2sigma or something similar)
    a : float
    alpha (distribution of relaxation times)
    """
    
    x = (1-a)*np.pi/2
    
    A = 2/(1-a)*np.arctanh(I*(1-a)*np.pi/np.arctanh(np.tan(x)))
    
    return None

def A1s(a):
    """
    Returns the value of the width parameter as fitted
    by D. Reta & N. Chilton. DOI: 10.1039/C9CP04301B
    
    Arguments
    ----------
    a : float
    alpha (distribution of relaxation times)
    
    Returns
    ----------
    A1s : float
    The value of the width parameter
    """

    A1s = 2.57*np.tan(alpha*np.pi/2)/alpha**0.096

    return A1s

def A2s(a):
    """
    Returns the value of the width parameter as fitted
    by D. Reta & N. Chilton. DOI: 10.1039/C9CP04301B
    
    Arguments
    ----------
    a : float
    alpha (distribution of relaxation times)
    
    Returns
    ----------
    A2s : float
    The value of the width parameter
    """

    A2s = 5.52*np.tan(alpha*np.pi/2)/alpha**0.30

    return A2s

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from scipy.integrate import quad
    
    f, ax = plt.subplots()
    s = np.linspace(-20,20,1001)
    
    alpha = 0.5
    
    I1s = 0.6827
    I2s = 0.9545
    
    A1s_own = A(I2s, alpha)
    A1s_ND = A1s(alpha)
    
    print(A1s_own, A1s_ND)
    
    #sigma=2.5
    #ax.plot(s, fuki_density(s, alpha), label='fuki')
    #ax.plot(s, 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-s**2/(2*sigma**2)), label='gauss')
    #ax.legend()
    #plt.show()