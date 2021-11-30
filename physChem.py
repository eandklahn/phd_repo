import numpy as np
#import scientificConstants as sc
import math
import scipy.constants as sc

kB = sc.Boltzmann

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

def qTm(m,p,T):
    
    mu = sc.physical_constants['atomic mass constant'][0]
    
    Vm = sc.R*T/p
    
    L = sc.h/np.sqrt(2*np.pi*m*mu*kB*T)
    
    qT = Vm/L**3
    return qT

def qR(B,s,T):
    """Returns the rotational partition function of a linear rotor
    Inputs
    B: rotational wavenumber(s) in cm^-1
    s: symmetry number
    T: temperature in K
    
    Outputs:
    qR: 
    """
    
    if isinstance(B, float):
        qR = kB*T/(s*sc.h*sc.c*10**2*B)
    elif isinstance(B, list):
        qR = 1/s*(kB*T/(sc.h*sc.c*10**2))**(3/2)*(np.pi/(B[0]*B[1]*B[2]))**(1/2)
    else:
        print('Type of B did not match any of the two cases')
        qR = None
    return qR

def qV(v,T):
    """Returns the vibrational partition function of a harmonic oscillator disregarding zero-point energy
    Inputs
    v: vibrational wavenumber in cm^-1
    T: temperature in K
    
    Outputs
    qV: vibrational partition function
    """
    
    kB = sc.Boltzmann
    
    qV = 1/(1-np.exp(-sc.h*sc.c*10**2*v/(kB*T)))
    
    return qV
    
def qV_tot(v_list, T):
    
    qV_returnval = 1
    for v in v_list:
        qV_returnval *= qV(v,T)
        
    return qV_returnval

def Uv(v,T):
    
    Na = sc.Avogadro
    kB = sc.Boltzmann
    
    Uv = Na*sc.h*sc.c*10**2*(v/(np.exp(sc.h*sc.c*10**2*v/(kB*T))-1))
    
    return Uv
    
def f_Ev(v,T):
    
    tv = sc.h*sc.c*10**2*v/kB
    f = (tv/T)**2*(np.exp(-tv/(2*T))/(1-np.exp(-tv/T)))**2
    
    return f
    
def mean_energy(Es, T):
    """
    Takes a Numpy array of energies "Es" and a temperature T
    and returns the mean energy based on Boltzmann populated states
    """
    
    exp_vector = np.exp(-sc.h*sc.c*10**2*Es/(sc.kB*T))
    E_mean = np.sum(Es*exp_vector)/np.sum(exp_vector)
    
    return E_mean.real
    
def U_v(v,T):
    
    return sc.N_A*(sc.h*sc.c*10**2*v)/(np.exp(sc.h*sc.c*10**2*v/(sc.Boltzmann*T))-1)
    