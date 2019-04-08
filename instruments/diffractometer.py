import numpy as np


class Diffractometer:

    def __init__(self):
        
        # Wavelength [Angstrom]
        self.L = 1
        # Temperature [Kelvin]
        self.T = 293
        # Magnetic field [Tesla]
        self.H = 1
        
        # Available angles
        self._range_omega = []
        self._range_chi = []
        self._range_phi = []
        self._range_2theta = []
    
    def _get_Qmax(self):
        
        _theta_max = np.radians(max(self._range_2theta)/2)
        
        Qmax = 4*np.pi*np.sin(_theta_max)/self.L
        
        return Qmax
    
    def _set_range_2theta(self,L):
    
        self._range_2theta = L
    
    def _observable_reflections(self, structure):
    
        Qmax = self._get_Qmax()
        
        h, k, l = 1, 1, 1
        while structure._Q_magnitude([h,0,0])<Qmax:
            h += 1
        while structure._Q_magnitude([0,k,0])<Qmax:
            k += 1
        while structure._Q_magnitude([0,0,l])<Qmax:
            l += 1
        hmax, kmax, lmax = h, k, l
        
        observables = []
        print('Calculating up to {} reflections'.format(2*hmax*2*kmax*2*lmax))
        for h in range(-hmax, hmax+1):
            for k in range(-kmax, kmax+1):
                for l in range(-lmax,lmax+1):
                    if (h,k,l) == (0,0,0):
                        continue
                    elif (h+k)%2!=0:
                        continue
                    elif structure._Q_magnitude([h,k,l])<Qmax:
                        observables.append((h,k,l))
            print(h)
        print(len(observables))
        
    
        