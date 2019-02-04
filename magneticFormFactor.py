import matplotlib.pyplot as plt
import numpy as np
from magnetism import landeG

def j_(n,s,D):
    """Calculates the integral <j_n(Q)> using the given analytical coefficients"""
    
    M = D['{}'.format(n)]['M']
    m = D['{}'.format(n)]['m']
    
    j = 0
    for t in range(len(m)):
        j += M[t]*np.exp(-m[t]*s**2)
    j += M[-1]
    if n>0:
        j *= s**2
    return j
    
def magF_(D, s, dipole=True):
    """Calculates the magnetic form factor for the compound with the coefficients given in D"""
    g = landeG(D['S'], D['L'], D['J'])

    f = 0
    if dipole:
        f = j_(0,s,D)+(2-g)/g*j_(2,s,D)
        return f
    else:
        return None

        
if __name__ == '__main__':
    # Magnetic form factor for DyIII as given in the input file for CHILSQ
    s = np.array([0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60])
    f = np.array([0.999900,0.988567,0.955746,0.904705,0.840101,0.767059,0.690328,0.539968,
                0.407135,0.298064,0.212983,0.149045,0.102148,0.068380,0.044657,0.028682,
                0.018692,0.013254,0.011191,0.011553])
                
    # Analytical approximation constants for DyIII
    DyIII = {
                '0': {'M':[0.1157, 0.3270, 0.5821, -0.0249], 'm': [15.0732, 6.7991,3.0202]},
                '2': {'M':[0.2523, 1.0914, 0.9345, 0.0250], 'm':[18.5172, 6.7362, 2.2082]},
                'J': 15/2, 'L': 5, 'S': 5/2
            }
            
    plt.plot(s,f)
    s = np.linspace(0,1.6,1000)
    plt.plot(s,magF_(DyIII, s))
    plt.xlabel(r'$\frac{\sin \theta}{\lambda}$ [$\AA^{-1}]$')
    plt.ylabel('Dysprosium mag. form factor')
    plt.show()