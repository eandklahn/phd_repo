import numpy as np
import scientificConstants as sc

print('Convert between Angstroms and keV\n')
print('Options:\n[1] Convert from keV to Angstroms\n[2] Convert from Angstroms to keV\n[3] Quit\n')

C = sc.h*sc.c/(10**-7*sc.e_charge)

while True:
    
    try:
        userIn = int(input())
        
        if userIn == 1:
            
            userIn = float(input('Insert an energy in keV: '))
            print('{:10.7f} keV corresponds to {:10.7f} Angstrom'.format(userIn, C/userIn))
            
        elif userIn == 2:
        
            userIn = float(input('Insert a wavelength in Angstrom: '))
            print('{:10.7f} Angstrom corresponds to {:10.7f} keV'.format(userIn, C/userIn))
            
        elif userIn == 3:
            
            break
            
    except:
        pass
            
    