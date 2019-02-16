import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

def _calculate_flipping_ratio(Fn, Fm_perp_Q, P, _e=1):
    """
    Calculates flipping ratios based on the equations found in ...
    
    Inputs
    Fn: nuclear structure factor for a reflection
    Fm_perp_Q: magnetic structure factor perpendicular to scattering vector in an orthonormal coordinate system
    P: polarization vector in same coordinates as Fm_perp_Q
    _e: flipper efficiency
    
    Outputs
    Tuple of (R_num R_den, R)
    R_num: numerator of the flipping ratio
    R_den: denominator of the flipping ratio
    R: R_num/R_den. The flipping ratio
    """
    
    Fn2 = Fn*np.conjugate(Fn)
    Fm2 = np.dot(np.transpose(Fm_perp_Q), np.conjugate(Fm_perp_Q))
    PFm = np.dot(np.transpose(P), Fm_perp_Q)
    
    R_num = Fn2+Fm2+2*(PFm*Fn)
    R_den = Fn2+Fm2-2*(PFm*Fn)*_e
    R = R_num/R_den
    
    return (float(R_num.real), float(R_den.real), float(R.real))
    
def _plot_peak(structure, h, k, l, _magnetic_field, _polarization):
    """Show the plot of a Bragg-peak specified by (h,k,l) using the nuclear and
    magnetic structure read from input structure.
    
    Inputs
    structure: an instance of crystallographyClasses.crystalStructure
    with a specified nuclear and magnetic structure
    hkl: the reciprocal lattice point for which to plot the peak
    _magnetic_field: strength of the magnetic field [unit: tesla]
    _polarization: number between 0 and 1 indicating beam polarization
    
    Outputs
    """
    
    # Nuclear structure independent of direction of magnetic field
    Fn = structure._calculate_structureFactor(h,k,l)
    
    t, p = np.linspace(0,np.pi,10), np.linspace(0,2*np.pi,20)
    T, P = np.meshgrid(t, p)
    
    X, Y, Z = np.sin(T)*np.cos(P), np.sin(T)*np.sin(P), np.cos(T)
    
    # Filling an array FR with magnetic field-dependent flipping ratios
    FR1 = np.zeros(P.shape)
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            t, p = T[i,j], P[i,j]
            Uijk = np.array([[np.sin(t)*np.cos(p)],[np.sin(t)*np.sin(p)],[np.cos(t)]])
            Fm_perp_Q = structure._calculate_Fm_perp_Q(h,k,l,-_magnetic_field*Uijk)
            FR1[i,j] = _calculate_flipping_ratio(Fn, Fm_perp_Q, _polarization*Uijk)[2]
    
    # Constructing figure from X, Y, Z and FR arrays
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    norm = matplotlib.colors.SymLogNorm(1, vmin=FR1.min(), vmax=FR1.max())
    
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                                       facecolors=cm.coolwarm(norm(FR1)),
                                       antialiased=False,
                                       shade=False)
                                        
    surf.set_edgecolors('None')
    m = cm.ScalarMappable(cmap=cm.coolwarm)
    m.set_array(FR1)
    fig.colorbar(m)
    
    # Calculating coordinates of unit cell axes in orthonormal system
    _a_coords_ijk = np.matmul(structure.IJK_Mct_ABC, np.array([[1],[0],[0]]))
    _b_coords_ijk = np.matmul(structure.IJK_Mct_ABC, np.array([[0],[1],[0]]))
    _c_coords_ijk = np.matmul(structure.IJK_Mct_ABC, np.array([[0],[0],[1]]))
    
    _a_coords_ijk = _a_coords_ijk/np.linalg.norm(_a_coords_ijk)*1.3
    _b_coords_ijk = _b_coords_ijk/np.linalg.norm(_b_coords_ijk)*1.3
    _c_coords_ijk = _c_coords_ijk/np.linalg.norm(_c_coords_ijk)*1.3
    
    # Setting text on figure indicating directions of crystal axes
    ax.text(1.2,0,0,'i')
    ax.text(0,1.2,0,'j')
    ax.text(0,0,1.2,'k')
    
    ax.text(_a_coords_ijk[0,0],_a_coords_ijk[1,0],_a_coords_ijk[2,0],'a')
    ax.text(_b_coords_ijk[0,0],_b_coords_ijk[1,0],_b_coords_ijk[2,0],'b')
    ax.text(_c_coords_ijk[0,0],_c_coords_ijk[1,0],_c_coords_ijk[2,0],'c')
    
    ax.set_aspect('equal')
    ax.set_axis_off()

    plt.show()
    
    
    
    