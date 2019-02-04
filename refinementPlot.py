from calculateCHILSQRs import *
import numpy as np
import matplotlib.pyplot as plt
from crystallography import *
from projectData import PND_susceptibility as D
import os

def plotLIScontent(ax, L, params, D, printBool=False):
    """
    ax: matplotlib.pyplot.Axes instance
    L: wavelength of the used radiation
    params: unit cell parameters of the crystal
    D: dictionary returned by calculateCHILSQRs.readCHILSQ_FRs
    """
    
    fileName = D['fileName']
    FRobs = D['FRobs']
    FRcalc = D['FRcalc']
    w = D['weight']
    q = D['q']
    s = D['s']
    nVar = D['nVar']
    
    R1 = calcR1(FRobs, FRcalc, doPrint=printBool)*100
    R2 = calcR2(FRobs, FRcalc, doPrint=printBool)*100
    R3 = calcR3(FRobs, FRcalc, w, doPrint=printBool)*100
    R4 = calcR4(FRobs, FRcalc, w, doPrint=printBool)*100
    X2 = calcXsqrd(FRobs, FRcalc, w, nVar)
    
    ax.text(0.01, 0.4, fileName+'\nR1 = {:3.2f}%\nR2 = {:3.2f}%\nR3 = {:3.2f}%\nR4 = {:3.2f}%\nX2 = {:3.2f}'.format(R1, R2, R3, R4, X2), transform=ax.transAxes, size=8)
    ax.errorbar(q, FRobs, yerr=s, c='k', marker='*', linestyle='-', label='Observed', linewidth=1)
    points, = ax.plot(q, FRcalc, 'r*-', label='Calculated', linewidth=1)
    ax.plot(q, FRcalc-FRobs, 'b-', label='Difference', linewidth=1)
    ax.legend(loc='upper right', fontsize=8)
    #ax.set_xlabel(r'$q = \frac{4\pi\sin{\theta}}{\lambda}$')
    #ax.set_ylabel(r'$FR = \frac{I_+}{I_-}$')
    
    return points,
    
if __name__ == '__main__':

    params = D.Dy_DBM_params
    L = 1.4
    
    fig = plt.figure()
    axes = []
    printing = False
    
    fileNames = [x for x in os.listdir() if x.endswith('.lis')]
    
    for n in range(len(fileNames)):
        ax1 = fig.add_subplot(len(fileNames), 1, n+1)
        axes.append(ax1)
        
        fileName = fileNames[n]
        data = readCHILSQ_FRs(fileNames[n])
        extFile = [x for x in list(data.keys()) if x.endswith('.ext')][0]
        A = data[extFile]
        
        rows_cols = A.shape
        h, k, l = A[:,0], A[:,1], A[:,2]
        FRobs = A[:,3]
        FRcalc = A[:,4]
        w = A[:,-1]
        q = q_from_theta(theta_from_d(d_from_hkl(h,k,l,params),L),L)
        
        sort_indices = q.argsort()
        q = q[sort_indices]
        FRobs = FRobs[sort_indices]
        FRcalc = FRcalc[sort_indices]
        w = w[sort_indices]
        R1 = calcR1(FRobs, FRcalc, doPrint=printing)*100
        R2 = calcR2(FRobs, FRcalc, doPrint=printing)*100
        R3 = calcR3(FRobs, FRcalc, w, doPrint=printing)*100
        R4 = calcR4(FRobs, FRcalc, w, doPrint=printing)*100
        X2 = calcXsqrd(FRobs, FRcalc, w, data['variables'])
        
        ax1.text(0.01, 0.4, extFile+'\nR1 = {:3.2f}%\nR2 = {:3.2f}%\nR3 = {:3.2f}%\nR4 = {:3.2f}%\nX2 = {:3.2f}'.format(R1, R2, R3, R4, X2), transform=ax1.transAxes, size=8)
        ax1.plot(q, FRobs, 'k*-', label='Observed', linewidth=1)
        ax1.plot(q, FRcalc, 'r*-', label='Calculated', linewidth=1)
        ax1.plot(q, FRcalc-FRobs, 'b-', label='Difference', linewidth=1)
        ax1.legend(loc='upper right', fontsize=8)
        ax1.set_xlabel(r'$q = \frac{4\pi\sin{\theta}}{\lambda}$')
        ax1.set_ylabel(r'$FR = \frac{I_+}{I_-}$')
        
        
    plt.show()




