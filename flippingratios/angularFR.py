import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import os
import time
from subprocess import Popen, PIPE
from calculateCHILSQRs import readCHILSQ_FRs

def plot_sphere_with_colors(H, WVLN, Pol, MAGF, ORIENTATION, thetaNum, phiNum):
    """
    Adapted from https://stackoverflow.com/questions/24218543/colouring-the-surface-of-a-sphere-with-a-set-of-scalar-values-in-matplotlib#
    """
    
    print('Calculating {} positions on sphere'.format(thetaNum*phiNum))
    
    theta = np.linspace(0, np.pi, num=thetaNum)
    phi   = np.linspace(0, 2*np.pi, num=phiNum)

    T, P = np.meshgrid(theta, phi)
    
    h, k, l = H
    
    fig = plt.figure()
    #ax1 = fig.add_subplot(1,2,1,projection='3d')
    ax2 = fig.gca(projection='3d')
    
    FR1 = simulate_FR_angle_dep([[h,k,l]], WVLN, ORIENTATION, Pol, MAGF, (thetaNum, phiNum))
    FR = FR1.ravel()
    
    X1, Y1, Z1 = np.sin(T)*np.cos(P), np.sin(T)*np.sin(P), np.cos(T)
    X, Y, Z = X1.ravel(), Y1.ravel(), Z1.ravel()
    
    norm = matplotlib.colors.SymLogNorm(1, vmin=FR.min(), vmax=FR.max())
    
    surf = ax2.plot_surface(X1, Y1, Z1, cmap=cm.coolwarm,
                                        facecolors=cm.coolwarm(norm(FR1)),
                                        antialiased=False,
                                        shade=False)
                                        
    surf.set_edgecolors('None')
    m = cm.ScalarMappable(cmap=cm.coolwarm)
    m.set_array(FR1)
    fig.colorbar(m)
    
    ax2.text(1.2,0,0,'i')
    ax2.text(0,1.2,0,'j')
    ax2.text(0,0,1.2,'k')
    
    ax2.set_aspect('equal')
    ax2.set_axis_off()
    
    plt.show()
    
def simulate_FR_angle_dep(HKL, WVLN, ORIENTATION, PVAL, MAGF, incrementsTuple, axes=None):

    numTheta = incrementsTuple[0]
    numPhi = incrementsTuple[1]
    
    Theta = np.linspace(0, np.pi, numTheta)
    Phi = np.linspace(0, 2*np.pi, numPhi)
    
    T, P = np.meshgrid(Theta, Phi)
    FR = np.zeros(P.shape)
    
    dim1, dim2 = P.shape
    tuplesToCheck = [(x,y) for x in range(dim1) for y in range(dim2)]
    
    niter, SF = 0, 0
    skippedFiles = []
    
    for set in tuplesToCheck:
            
        t, p = T[set[0],set[1]], P[set[0],set[1]]
        
        rflname = '{}.ext'.format(niter)
        
        magFieldDirecion = [np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)]

        POLAR = magFieldDirecion + [-PVAL, PVAL]
        
        make_rfl(WVLN, ORIENTATION, POLAR, MAGF, HKL, NAME=rflname)
        
        while True:
            try:
                PE = runchilsq('DBM.cry', rflname)
                FRcalc = readCHILSQ_FRs('chilsq.lis')
                FRcalc = FRcalc[rflname]['FRcalc'][0]
                break
            except KeyError:
                continue
        
        FR[set[0], set[1]] = FRcalc
            
        try:
            os.remove(rflname)
        except PermissionError:
            skippedFiles.append(niter)
        
        niter += 1
        if PE: SF += 1
        percentageDone = niter/(dim1*dim2)*100
        
        print('PROGRESS: |'+'#'*int(percentageDone)+' '*(100-int(percentageDone))+
              '| {:>6.2f}% niter={}, skipped files={}'.format(round(percentageDone,2), niter, SF), end='\r')
        
    print('\nRemoving a total of {} skipped file(s)'.format(len(skippedFiles)))
    for file in skippedFiles:
        os.remove('{}.ext'.format(file))
        
    if axes is not None:
        
        axes.plot([0, 0], [np.pi/2, np.pi/2], [0.5, 1.5], c='k') # Plotting x-line
        axes.text(0,np.pi/2,1.6,'i')
        axes.plot([np.pi/2, np.pi/2], [np.pi/2, np.pi/2], [0.5, 1.5], c='k') # Plotting x-line
        axes.text(np.pi/2,np.pi/2,1.6,'j')
        axes.plot([0, 0], [0, 0], [0.5, 1.5], c='k') # Plotting x-line
        axes.text(0,0,1.6,'k')
        
        surf = axes.plot_surface(P, T, FR, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
        
        axes.set_title('FR for ({},{},{})'.format(HKL[0][0], HKL[0][1], HKL[0][2]))
        axes.set_xlabel(r'$\phi$')
        axes.set_ylabel(r'$\theta$')
        axes.set_zlabel(r'$FR$')
        axes.set_zlim(0.5,1.5)
        
        axes.set_aspect('equal')
    
    return FR

def runchilsq(cryFile, file):
    
    PE = False
    
    remove_first = ['chilsq.lis']
    for item in remove_first:
        while True:
            try:
                os.remove(item)
                break
            except PermissionError:
                PE = True
                time.sleep(0.015)
                print('Caught PermissionError')
                continue
            except FileNotFoundError:
                break

    string = '{}\n{} 5 5 1\n\n\n'.format(cryFile, file)

    chilsq = Popen('chilsq', stdin=PIPE, stdout=PIPE, stderr=PIPE)
    chilsq.stdin.write(bytes(string, 'utf8'))
    output = chilsq.communicate()
    
    return PE

def make_rfl(WVLN, ORIENTATION, POLAR, MAGF, HKL, NAME=None):
    
    """
    A good explanation of the content of the data file was found here: https://www.ill.eu/sites/ccsl/mk4man/c5node5.html
    Most important information given in the file (in the case the page disappears) is
    
    Polarisation followed by 5 numbers: the Up and Down polarising efficiencies and the polarisation direction in
    orthogonal crystal coordinates.
    
    Wavelength followed by its value and, if needed, the absorption coefficient.
    
    Orientation followed by 9 components of the matrix giving the directions of the diffractometer axes in orthogonal
    crystal coordinates.
    
    NOTE an important difference: In CHILSQ, the first 3 numbers in polarisation are coordinates. After that comes
    polarisation values.
    """
    
    if NAME is not None:
        f = open(NAME, 'w')
    else:
        f = open('rfl.ext', 'w')
    
    # First line in the ext-file is the wavelength of radiation in Angstrom
    f.write('#Wavelength{:>10.4f}\n'.format(
    WVLN))
    
    # Second line is the orientation matrix (U) which gives the direction of the diffractometer axes in orthogonal crystal
    # coordinates. Probably in CCSL, the IJK-system is used as the orthogonal crystal coordinate system.
    f.write('#Orientation{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(
    *ORIENTATION))
    
    # The polarisation is given in orthogonal crystal coordinates. Again, the orthogonal crystal coordinates used here
    # are probably the coordinates in the IJK-system.
    f.write('#Polarisation{:>10.4f}{:>10.4f}{:>10.4f}{:>10.4f}{:>10.4f}\n'.format(
    *POLAR))
    
    # Fourth line is the strength of the magnetic field.
    f.write('#Magnetic_field  {:<5.4f}\n'.format(
    MAGF))
    
    for set in HKL:
        f.write('{:>5d}{:>5d}{:>5d}{:>10.4f}{:>10.4f}{:>10.4f}{:>10.4f}{:>10.4f}{:>10.4f}{:>5d}\n'.format(
        set[0], set[1], set[2], 1, 0.01, 1, 1, 1, 1, 1))
    
    f.close()

def calculate_FR(tuple):

    WVLN        = tuple[0]
    ORIENTATION = tuple[1]
    POLARVAL    = tuple[2]
    MAGF        = tuple[3]
    
    t, p        = tuple[4], tuple[5]
    h, k, l     = tuple[6], tuple[7], tuple[8]
    
    x, y, z = np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)
    POLAR   = [x, y, z, -POLARVAL, POLARVAL]
    HKL     = [[h,k,l]]
    saveas  = '{}.ext'.format(int(int((t+p)*100)))
    
    make_rfl(WVLN, ORIENTATION, POLAR, MAGF, HKL, NAME=saveas)
    
    runchilsq('DBM.cry', saveas)
    
    FR = readCHILSQ_FRs('chilsq.lis')
    FR = FR[saveas]['FRcalc'][0]
    os.remove(saveas)
    
    f = open(saveas[:-4]+'.res', 'w')
    f.write('{} {} {}\n'.format(t, p, FR))
    f.close()
    
    
    
if __name__ == '__main__':
    
    WVLN = 1.4
    
    ORIENTATION = [-0.2056,  0.8439,  0.4956, -0.7930, -0.4404,  0.4208,  0.5734, -0.3065,  0.7598]
    
    P = 0.78
    POLAR = [0.4956,    0.4208,    0.7598, -P, P]
    
    MAGF = 0.5
    
    thetaNum = 10
    phiNum = 2*thetaNum
    
    """This is where it is actually possible to run a sphere plotting"""
    plot_sphere_with_colors([4,4,4], WVLN, P, MAGF, ORIENTATION, thetaNum, phiNum)
    
    """Below is a start at an attempt to parallelise the process. calculate_FR is the important function together
    with listToGive which is a list comprehension to give a list of tuples that can be iterated through."""
    #h, k, l = 1, 1, 1
    #t = np.linspace(0,np.pi,thetaNum)
    #p = np.linspace(0,2*np.pi,phiNum)
    #
    #listToGive = [(WVLN, ORIENTATION, P, MAGF, theta, phi, h, k, l) for theta in t for phi in p]
    #for item in listToGive:
    #    calculate_FR(item)
    