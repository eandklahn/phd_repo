import numpy as np
import matplotlib.pyplot as plt
from cifoperations import mj0, mj2, mj4, mj6

def magF_(ion, s, L=None, S=None, J=None, type='j0'):
    """
    Calculates the magnetic form factor of 'ion'
    s: sinT/lambda
    """
    
    j0 = _calculate_exponential_approx(ion, s, type='j0')
    
    f = 0
    if type == 'j0':
        return j0
    elif type == 'dipole':
        g = landeG(L,S,J)
        j2 = _calculate_exponential_approx(ion, s, type='j2')
        # Use expression from McPhase manual for magnetic form factor
        f = j0+(2-g)/g*j2
        return f
    else:
        return None

def landeG(L,S,J):
    """Function that returns the Landé g-factor for quantum numbers L, S and J"""
    
    g = 3/2 + (S*(S+1)-L*(L+1))/(2*J*(J+1))
    
    return g
    
def _plot_formfactor(ion, type, Sstart, Send):

    s = np.linspace(Sstart, Send, 1000)
    f = _calculate_exponential_approx(ion, s, type=type)
    
    plt.plot(s,f)
    plt.show()
    
def _calculate_exponential_approx(ion, s, type='elec'):
        
        dict = {}
        if type == 'elec':
            dict = factabe.elec
            n = 0
        elif type == 'j0':
            dict = mj0.mj0
            n = 0
        elif type == 'j2':
            dict = mj2.mj2
            n = 2
        
        if ion in dict.keys():
            
            a = dict[ion]['a']
            b = dict[ion]['b']
            c = dict[ion]['c']
            
            f = 0
            for i in range(len(a)):
                f += a[i]*np.exp(-b[i]*s**2)
            f += c
            if n>1:
                f *= s**2
            
            return f
            
        else:
            print('The requested ion is not in the form factor table')

def _update_magj0formfactortable():

    f = open('_j0_ form factors for rare earth ions.html', 'r')
    data = f.readlines()
    f.close()
    
    data = data[47:244]
    
    ions = []
    
    n=0
    while n<len(data):
        if '</tr>' in data[n]:
            n += 1
        elif 'LEFT' in data[n]:
            line = data[n]
            ion = line.split('>')[2].split('<')[0]
            iondata = [ion]
            for i in range(1,8):
                D = data[n+i].split('>')[1].split('<')[0]
                iondata.append(D)
            ions.append(iondata)
            n += 8
        else:
            n += 1
    
    f = open('magformfactorj0table.py', 'w')
    
    f.write('''"""
These form factors have been obtained from
https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html
"""''')
    
    f.write('\n\n')
    
    f.write('magj0 = {\n')
    
    for iondata in ions:
        
        ion = iondata[0]
        a = [float(x) for x in iondata[1:-1:2]]
        b = [float(x) for x in iondata[2:-1:2]]
        c = float(iondata[-1])
        
        texta = "'a': [{:<9.6},{:<9.6},{:<9.6}]".format(a[0], a[1], a[2])
        textb = "'b': [{:<9.6},{:<9.6},{:<9.6}]".format(b[0], b[1], b[2])
        textc = "'c':  {}".format(c)
        
        f.write("'{}':".format(ion)+' '*(10-len(ion)-3)+'{'+texta+',\n           '+textb+',\n           '+textc+'},\n\n')
        
    f.write('}')
    f.close()
    
    print('Magnetic form factor table j0 updated')

def _update_magj2formfactortable():

    f = open('_j2_ form factors for rare earth atoms and ions.html', 'r')
    data = f.readlines()
    f.close()
    
    data = data[47:244]
    
    ions = []
    
    n=0
    while n<len(data):
        if '</tr>' in data[n]:
            n += 1
        elif 'LEFT' in data[n]:
            line = data[n]
            ion = line.split('>')[2].split('<')[0]
            iondata = [ion]
            for i in range(1,8):
                D = data[n+i].split('>')[1].split('<')[0]
                iondata.append(D)
            ions.append(iondata)
            n += 8
        else:
            n += 1
    
    f = open('magformfactorj2table.py', 'w')
    
    f.write('''"""
These form factors have been obtained from
https://www.ill.eu/sites/ccsl/ffacts/ffachtml.html
"""''')
    
    f.write('\n\n')
    
    f.write('magj2 = {\n')
    
    for iondata in ions:
        
        ion = iondata[0]
        a = [float(x) for x in iondata[1:-1:2]]
        b = [float(x) for x in iondata[2:-1:2]]
        c = float(iondata[-1])
        
        texta = "'a': [{:<9.6},{:<9.6},{:<9.6}]".format(a[0], a[1], a[2])
        textb = "'b': [{:<9.6},{:<9.6},{:<9.6}]".format(b[0], b[1], b[2])
        textc = "'c':  {}".format(c)
        
        f.write("'{}':".format(ion)+' '*(10-len(ion)-3)+'{'+texta+',\n           '+textb+',\n           '+textc+'},\n\n')
        
    f.write('}')
    f.close()
    
    print('Magnetic form factor table j2 updated')
    
def _update_elecformfactortable():

    filename = 'Atomic form factors.html'
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    f = open('elecformfactortable.py', 'w')
    
    f.write('''"""
These form factors have been obtained from
http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
which uses the values calculated in http://it.iucr.org/Cb/ch6o1v0001/
"""''')
    
    f.write('\n\nelec = {\n')
    
    formfactors = data[358:570]
    
    for line in formfactors:
        line = line.split(';')
        line = [x.split('=') for x in line]
        
        ion = line[0][1][1:-1]
        a1 = float(line[1][1])
        b1 = float(line[2][1])
        a2 = float(line[3][1])
        b2 = float(line[4][1])
        a3 = float(line[5][1])
        b3 = float(line[6][1])
        a4 = float(line[7][1])
        b4 = float(line[8][1])
        c  = float(line[9][1])
        
        
        
        texta = "'a': [{:<9.6},{:<9.6},{:<9.6},{:<9.6}]".format(a1, a2, a3, a4)
        textb = "'b': [{:<9.6},{:<9.6},{:<9.6},{:<9.6}]".format(b1, b2, b3, b4)
        textc = "'c':  {}".format(c)
        
        f.write("'{}':".format(ion)+' '*(10-len(ion)-3)+'{'+texta+',\n           '+textb+',\n           '+textc+'},\n\n')
    
    f.write('}')
        
    f.close()
    
    print('Electronic form factor table updated')
    
if __name__ == '__main__':
    
    s = np.linspace(0,10,101)
    magF_(ion, s, L=5, S=2/5, J=15/2, type='dipole')
    
    
    #_update_elecformfactortable()
    #_update_magj0formfactortable()
    #_update_magj2formfactortable()
    #
    #_plot_formfactor('Co2')