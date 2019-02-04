import numpy as np
import scientificConstants as sc
import math
"""Functions for the calculation of exercises in Physical Chemistry I"""

def psi_box(x, n, L):
    """Wave function for the box potential with infinite walls"""
    return np.sqrt(2/L)*np.sin(n*np.pi*x/L)
    
def hermitePolynomial(n,x):
    """Returns the value of the n'th Hermite polynomial according to the recursion H_(n+1) = 2x H_n - 2nH_(n-1)"""
    
    assert n>-1
    
    if n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return 2*x*hermitePolynomial(n-1,x) - 2*(n-1)*hermitePolynomial(n-2,x)
        
def psi_harmonic(x,n,m,k):
    """
    Returns the value of the normalized n'th vibrational wavefunction in the harmonic oscillator
    x: position(s) at which the function is evaluated
    n: vibrational quantum number
    m: reduced mass of the vibration (in atomic mass units)
    k: force constant for the vibration (in N/m)
    """
    
    m *= sc.amu
    w = np.sqrt(k/m)
    epsilon = np.sqrt(m*w/sc.hBar)*x
    
    return (m*w/(np.pi*sc.hBar))**(1/4)*1/np.sqrt(2**n*math.factorial(n))*hermitePolynomial(n,epsilon)*np.exp(-epsilon**2/2)
    
def qV(v,T):
    """Returns the vibrational partition function of a harmonic oscillator disregarding zero-point energy
    Inputs
    v: vibrational wavenumber in cm^-1
    T: temperature in K
    
    Outputs
    qV: vibrational partition function
    """
    
    qV = 1/(1-np.exp(-sc.h*sc.c*10**2*v/(sc.kB*T)))
    
    return qV
    
def qV_tot(v_list, T):

    qV_returnval = 1
    for v in v_list:
        qV_returnval *= qV(v,T)
        
    return qV_returnval

def Uv(v,T):

    Uv = sc.Na*sc.h*sc.c*10**2*(v/(np.exp(sc.h*sc.c*10**2*v/(sc.kB*T))))
    
    return Uv
    
def f_Ev(v,T):

    tv = sc.h*sc.c*10**2*v/sc.kB

    f = (tv/T)**2*(np.exp(-tv/(2*T))/(1-np.exp(-tv/T)))**2
    
    return f