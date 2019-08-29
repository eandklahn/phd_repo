import numpy as np
from numpy import cos, sin, arcsin, radians, sqrt, pi as PI
from cifoperations import crystallographyClasses as cc

def rota(a):

    a = np.radians(a)

    R = np.mat([[1, 0,         0         ],
                [0, np.cos(a), -np.sin(a)],
                [0, np.sin(a), np.cos(a) ]])
    
    return R

def rotb(b):

    b = np.radians(b)

    R = np.mat([[np.cos(b), 0, -np.sin(b)],
                [0,         1, 0         ],
                [np.sin(b), 0, np.cos(b) ]])
    
    return R
    
def rotc(c):

    c = np.radians(c)

    R = np.mat([[np.cos(c), -np.sin(c), 0],
                [np.sin(c), np.cos(c),  0],
                [0,         0,          1]])
    
    return R
    
def d_from_hkl(h,k,l,params):
    """Calculates the d-spacing corresponding to a given set of hkl-values based on the formula for the triclinic lattice in Giacovazzo, table 2.2
    Inputs
    h,k,l: reciprocal lattice point
    params: a list of crystallographic lattice constants [a,b,c,alpha,beta,gamma]. Angles must be in degrees.
    
    Outputs
    d: d-spacing corresponding to a given reciprocal lattice point
    """
    
    # Extracting lattice constants from params
    a = params[0]
    b = params[1]
    c = params[2]
    A = radians(params[3])
    B = radians(params[4])
    Y = radians(params[5])
    
    one_over_d_dhkl_squared = (1-cos(A)**2-cos(B)**2-cos(Y)**2 + 2*cos(A)*cos(B)*cos(Y))**(-1)*((h/a*sin(A))**2 + (k/b*sin(B))**2 + (l/c*sin(Y))**2 \
                               + 2*k*l/(b*c)*(cos(B)*cos(Y) - cos(A)) + 2*l*h/(c*a)*(cos(Y)*cos(A) - cos(B)) + 2*h*k/(a*b)*(cos(A)*cos(B) - cos(Y)))
                               
    d = sqrt(1/one_over_d_dhkl_squared)
    
    return d
    
def theta_from_d(d, L):
    """Calculates theta using Braggs law
    Inputs
    d: lattice spacing for which the scattering angle should be calculated
    L: wavelength of the used radiation
    
    Outputs
    t: theta in radians"""
    
    t = arcsin(L/(2*d))
    
    return t
    
def q_from_theta(t, L):
    """Calculates the magnitude of the scattering vector q using q=4*pi*sin(theta)/lambda
    Inputs
    t: scattering angle in radians
    L: wavelength of used radiation
    
    Outputs
    q: magnitude of scattering vector
    """
    
    q = 4*PI*sin(t)/L
    
    return q

def q_from_hkl(h,k,l,L,params):
    """Combines previous functions to calculate a q-value directly from a**2
    given set of hkl-values and a wavelength
    Inputs
    h,k,l: reciprocal lattice point
    L: wavelength of used radiation in Angstrom
    params: a list of crystallographic lattice constants
    
    Outputs
    q: magnitude of scattering vector
    """
    
    d = d_from_hkl(h,k,l,params)
    t = theta_from_d(d,L)
    q = q_from_theta(t,L)
    
    return q