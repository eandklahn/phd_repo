import matplotlib.pyplot as plt
from calculateCHILSQRs import readCHILSQ_FRs
import numpy as np

def compare_both():

    D1 = readCHILSQ_FRs('chilsq_own.lis')
    D2 = readCHILSQ_FRs('chilsq_old.lis')
    
    f, ax = plt.subplots(nrows=(len(list(D1.keys()))-1), ncols=2)
    ax[0,0].errorbar(x=list(range(len(D1['out73.ext']['FRobs']))),
                    y=D1['out73.ext']['FRobs'],
                    fmt='k*',
                    yerr=np.sqrt(2/D1['out73.ext']['weight']),
                    label='measured')
    ax[0,0].plot(D1['out73.ext']['FRcalc'], 'r*', label='calculated')
    ax[1,0].errorbar(x=list(range(len(D1['out113.ext']['FRobs']))),
                    y=D1['out113.ext']['FRobs'],
                    fmt='k*',
                    yerr=np.sqrt(2/D1['out113.ext']['weight']),
                    label='measured')
    ax[1,0].plot(D1['out113.ext']['FRcalc'], 'r*', label='calculated')
    ax[0,1].errorbar(x=list(range(len(D2['rfl1.ext']['FRobs']))),
                    y=D2['rfl1.ext']['FRobs'],
                    fmt='k*',
                    yerr=np.sqrt(2/D2['rfl1.ext']['weight']),
                    label='measured')
    ax[0,1].plot(D2['rfl1.ext']['FRcalc'], 'r*', label='calculated')
    ax[1,1].errorbar(x=list(range(len(D2['rfl2.ext']['FRobs']))),
                    y=D2['rfl2.ext']['FRobs'],
                    fmt='k*',
                    yerr=np.sqrt(2/D2['rfl2.ext']['weight']),
                    label='measured')
    ax[1,1].plot(D2['rfl2.ext']['FRcalc'], 'r*', label='calculated')
    
    chi_sqrd = [1.43, 4.79, 4.48, 8.68]
    titles = ['orientation 1, new', 'orientation 1, old','orientation 2, new','orientation 2, old']
    for i,a in enumerate(f.get_axes()):
        a.set_title(titles[i] + ' ($\chi^2$ = {})'.format(chi_sqrd[i]))
        a.set_ylabel('flipping ratio')
        a.set_xlabel('flipping ratio number')
    ax[1,0].legend(loc='upper right')
    f.tight_layout()
    plt.show()
    
def show_one_result():

    D = readCHILSQ_FRs('chilsq_own.lis')
    
    O1_name = 'out73.ext'
    O2_name = 'out113.ext'
    
    O1 = D[O1_name]
    O2 = D[O2_name]
    
    O1_FRcalc = O1['FRcalc']
    O1_FRobs = O1['FRobs']
    O1_err = np.sqrt(2/O1['weight'])
    O1_x = np.arange(len(O1_FRcalc))
    
    O2_FRcalc = O2['FRcalc']
    O2_FRobs = O2['FRobs']
    O2_err = np.sqrt(2/O2['weight'])
    O2_x = np.arange(len(O2_FRcalc))
    
    # Plotting
    obs_fmt = 'k--*'
    calc_fmt = 'r--*'
    
    f, ax = plt.subplots(nrows=2, sharex=True)
    ax[0].errorbar(O1_x, O1_FRobs, yerr=O1_err, fmt=obs_fmt)
    ax[0].plot(O1_x, O1_FRcalc, calc_fmt)
    
    ax[1].errorbar(O2_x, O2_FRobs, yerr=O2_err, fmt=obs_fmt)
    ax[1].plot(O2_x, O2_FRcalc, calc_fmt)
    
    for a in f.get_axes():
        a.set_ylabel('FR = \frac{I_+}{I_-}')
        a.set_xlabel('Measurement number')
    
    plt.show()
    
if __name__ == '__main__':

    show_one_result()
    
    
    
    
    