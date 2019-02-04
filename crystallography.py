import numpy as np
from numpy import cos, sin, arcsin, radians, sqrt, pi as PI

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
    
class crystalStructure:
    """A class for making calculations and holding transformation matrices relating to individual crystal structures."""

    def __init__(self, params, atoms=None):

        self.cellParams = params
        self.atoms = []
        
        self.a =        self.cellParams[0]
        self.b =        self.cellParams[1]
        self.c =        self.cellParams[2]
        self.alpha =    np.radians(self.cellParams[3])
        self.beta =     np.radians(self.cellParams[4])
        self.gamma =    np.radians(self.cellParams[5])
        
        # Metric matrix for the compound
        self.G = np.mat([[self.a**2, self.a*self.b*np.cos(self.gamma), self.a*self.c*np.cos(self.beta)],
                         [self.a*self.b*np.cos(self.gamma), self.b**2, self.b*self.c*np.cos(self.alpha)],
                         [self.a*self.c*np.cos(self.beta), self.b*self.c*np.cos(self.alpha), self.c**2]])
        
        # Volume of the unit cell
        self.V = np.sqrt(np.linalg.det(self.G))
        
        self.aStar =        self.b*self.c*np.sin(self.alpha)/self.V
        self.bStar =        self.a*self.c*np.sin(self.beta)/self.V
        self.cStar =        self.a*self.b*np.sin(self.gamma)/self.V      
        self.alphaStar =    (np.cos(self.beta)*np.cos(self.gamma) - np.cos(self.alpha))/(np.sin(self.beta)*np.sin(self.gamma))
        self.betaStar =     (np.cos(self.alpha)*np.cos(self.gamma) - np.cos(self.beta))/(np.sin(self.alpha)*np.sin(self.gamma))
        self.gammaStar =    (np.cos(self.alpha)*np.cos(self.beta) - np.cos(self.gamma))/(np.sin(self.alpha)*np.sin(self.beta))
        
        # Volume of the reciprocal cell
        self.VStar = 1/self.V
        
        # Metric matrix of the reciprocal lattice
        self.GStar = self.G.I
        
        # r, q and p as they are defined in the text in my bachelor project
        r = np.cos(self.beta)
        
        q = ((np.cos(self.gamma)-np.cos(self.beta)*np.cos(self.alpha))/(np.sin(self.alpha)))
        
        p = np.sqrt(1-(q**2+r**2))
        
        # Matrix for basis transformation from orthonormal CCSL basis to crystallographic basis
        # This is the basis transformation matrix (called 'M' in Giacovazzo) from the basis ijk to the basis abc
        self.ABC_Mbt_IJK = np.mat([[self.a*p, self.a*q,                  self.a*r                 ],
                                   [0,        self.b*np.sin(self.alpha), self.b*np.cos(self.alpha)],
                                   [0,        0,                         self.c                   ]])
        
        # Matrix for coordinate transformation from orthonormal CCSL basis to crystallographic basis
        self.ABC_Mct_IJK = np.linalg.inv(np.transpose(self.ABC_Mbt_IJK))
        
        # Matrix for basis transformation from crystallographic basis to orthonormal CCSL basis
        self.IJK_Mbt_ABC = np.linalg.inv(self.ABC_Mbt_IJK)
        
        # Matrix for coordinate transformation from crystallographic basis to orthonormal CCSL basis
        self.IJK_Mct_ABC = np.linalg.inv(np.transpose(self.IJK_Mbt_ABC))
        
    def bondLength(self, A1, A2):
        """Calculates the length of a bond in the structure based on the metric matrix of the structure"""
    
        v = np.mat([[A1.position[0,0]-A2.position[0,0]],
                    [A1.position[1,0]-A2.position[1,0]],
                    [A1.position[2,0]-A2.position[2,0]]])
        
        L = v.T*self.metricMatrix*v
        
        return L

    def bondAngle(self, A1, A2, A3):
        """Calculates the angle between the bonds originating in A2 and going to A1 and A3."""
        pass
    
    def addAtom(self, A):
        """Adds the atom A to the atoms contained in the structure."""
        
        self.atoms.append(A)
        
    def removeAtom(self, label):
        """Removes the atom with the given label from the atom container of the structure."""
        pass
        
    def electronicStructureFactor(h,k,l):
        """Calculates the structure factor for a particular point in reciprocal space"""
        Q = h*self.alphaStar + k*self.betaStar + l*self.gammaStar
        F = 0
        for atom in atoms:
            x = atom.position[0]
            y = atom.position[1]
            z = atom.position[2]
            F += atom.formfactor(Q)*np.exp(-2*np.pi*j(h*x+k*y+l*z))
            
class Atom:

    def __init__(self, position, type, label=None, ionisation=None):
    
        self.position = np.mat(position)
        self.label = label
        
    def createLabel(self, structure):
        """Makes an appropriate label for the atom based on the other atoms in the structure"""
        pass