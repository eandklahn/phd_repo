import numpy as np
from numpy import cos, sin, arcsin, radians, sqrt, pi as PI
from cifoperations import cc

"""
Giacovazzo:
Giacovazzo, C., Fundamentals of crystallography. 3. ed. ed.; Oxford University Press: Oxford, 2011
"""

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

def calculate_reciprocal_lattice(a,b,c,A,B,C,tol=1e-10):
    """
    Calculates reciprocal lattice parameters based
    on the given direct lattice parameters
    
    a,b,c : floats
            unit cell side lengths
    A,B,C : floats
            alpha, beta, gamma angles in degrees
    
    Returns
    ----------
    t : tuple
        A tuple containing (a_, b_, c_, A_, B_, C_) which are
        a, b, c, alpha, beta, gamma for the reciprocal lattice
        with the angles in degrees
    """
    
    G = calculate_G_matrix(a,b,c,A,B,C)
    G_ = np.linalg.inv(G)
    
    a_ = np.sqrt(G_[0,0])
    b_ = np.sqrt(G_[1,1])
    c_ = np.sqrt(G_[2,2])
    
    cosA_ = G_[1,2]/(b_*c_)
    A_ = np.degrees(np.arccos(cosA_))
    
    cosB_ = G_[0,2]/(a_*c_)
    B_ = np.degrees(np.arccos(cosB_))
    
    cosC_ = G_[0,1]/(a_*b_)
    C_ = np.degrees(np.arccos(cosC_))
    
    # Old approach
    #V = np.sqrt(np.linalg.det(G))
    #a_ = b*c*np.sin(np.radians(A))/V
    #b_ = a*c*np.sin(np.radians(B))/V
    #c_ = a*b*np.sin(np.radians(C))/V
    #
    #A_ = np.degrees(np.arcsin(V/(a*b*c*np.sin(np.radians(B))*np.sin(np.radians(C)))))
    #B_ = np.degrees(np.arcsin(V/(a*b*c*np.sin(np.radians(A))*np.sin(np.radians(C)))))
    #C_ = np.degrees(np.arcsin(V/(a*b*c*np.sin(np.radians(A))*np.sin(np.radians(B)))))
    
    return (a_, b_, c_, A_, B_, C_)
    
def calculate_G_matrix(a,b,c,A,B,C):
    """
    Calculates the G-matrix (metric matrix) in Giacovazzo
    based on the given direct unit cell parameters
    
    a,b,c : floats
            unit cell side lengths
    A,B,C : floats
            alpha, beta, gamma angles in degrees
    
    Returns
    ----------
    G : array
        A 3x3-array with the G-matrix
    """
    
    G11 = a**2
    G12 = a*b*np.cos(np.radians(C))
    G13 = a*c*np.cos(np.radians(B))
    G22 = b**2
    G23 = b*c*np.cos(np.radians(A))
    G33 = c**2
    
    G = np.array([[G11, G12, G13],
                  [G12, G22, G23],
                  [G13, G23, G33]])
    
    return G

def calculate_B_matrix(a,b,c,A,B,C):
    """
    Calculates the B-matrix as given in Busing & Levy
    (Acta Cryst. (1967). 22, 457)
    based on the given direct unit cell parameters
    Reciprocal lattice parameters are given with _ after. (a_, b_, etc.)
    
    a,b,c : floats
            unit cell side lengths
    A,B,C : floats
            alpha, beta, gamma angles in degrees
    
    Returns
    ----------
    B : array
        A 3x3-array with the B-matrix
    """
    
    a_, b_, c_, A_, B_, C_ = calculate_reciprocal_lattice(a,b,c,A,B,C)
    
    B11 = a_
    B12 = b_*np.cos(np.radians(C_))
    B13 = c_*np.cos(np.radians(B_))
    B22 = b_*np.sin(np.radians(C_))
    B23 = -c_*np.sin(np.radians(B_))*np.cos(np.radians(A))
    B33 = 1/c
    
    B = np.array([[B11, B12, B13],
                  [  0, B22, B23],
                  [  0,   0, B33]])
    
    return B

def calculate_Euler_matrix(t1,t2,t3):
    """
    Calculates the coordinate transformation matrix associated with
    the Euler angles t1, t2 and t3.
    Reference: Giacovazzo, p. 78
    
    t1, t2, t3 : floats
                 angles in degrees
    
    Returns
    ----------
    R : array
        A 3x3-array corresponding to the transformation matrix
    """
    
    t1, t2, t3 = np.radians(t1), np.radians(t2), np.radians(t3)
    
    c = np.cos
    s = np.sin

    R_Eu = np.array([[ c(t1)*c(t3)-s(t1)*c(t2)*s(t3),  s(t1)*c(t3)+c(t1)*c(t2)*s(t3), s(t2)*s(t3)],
                     [-s(t1)*c(t2)*c(t3)-c(t1)*s(t3), -s(t1)*s(t3)+c(t1)*c(t2)*c(t3), s(t2)*c(t3)],
                     [                   s(t1)*s(t2),                   -c(t1)*s(t2),       c(t2)]])
                     
    return R_Eu
    
def transform_property(M, Q):
    """
    Transforms the property according to the rule in Giacovazzo,
    that Q'=M*Q*M.T, when M is the basis transformation from
    A to A'.
    
    M : array
        Basis transformation matrix from A to A'
        
    Q : array
        Property given in the A-basis
    
    Returns
    ----------
    Q_transformed : array
                    A 3x3-array giving the property Q in the A'-basis
    """
    
    Q_transformed = np.matmul(M, np.matmul(Q, np.transpose(M)))
    
    return Q_transformed

def dot_product(v1, v2, G=np.identity(3)):
    
    return float(np.matmul(v1.T, np.matmul(G, v2)))

def vector_norm(v, G=np.identity(3)):
    
    norm_sqrd = dot_product(v, v, G)
    
    return float(np.sqrt(norm_sqrd))
  
def angle_between(v1, v2, G=np.identity(3)):
    
    cosangle = dot_product(v1, v2, G)/(vector_norm(v1, G)*vector_norm(v2, G))
    angle = np.arccos(cosangle)
    angle = np.degrees(angle)
    
    return float(angle)

























