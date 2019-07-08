import numpy as np
import crystallography as c

def readCHILSQ_FRs(fileName, params=None):
    """
    Reads the CHILSQ output for the lists of observed and calculated structure factors.
    
    Input
    fileName: name of the file to be read (including extension)
    
    Output
    dictOut: a dictionary with reflection file names as keys and arrays of flipping ratio data as printed by CHILSQ.
             An additional key 'variables' gives the number of variables that was refined during that run of the program to be used in calculating chi-squared values.
    """

    f = open(fileName, 'r')
    D = f.readlines()
    f.close()
    
    D = [line.strip() for line in D]
    N = len(D)
    
    names = []
    beginnings = []
    endings = []
    data = []
    variables = 0
    
    n = 0
    while n < N:
        line = D[n]
        if 'Input file' in line:
            names.append(line.split()[2])
        if 'h    k    l' in line:
            beginnings.append(n+1)
            while len(D[n]) > 2:
                n += 1
            endings.append(n)
        if 'basic variable(s)' in line:
            variables = float(line.split()[0])
        n += 1
    
    print(beginnings)
    
    dictOut = {}
    for n in range(len(names)):
        A = np.array([line.split() for line in D[beginnings[n]:endings[n]]], dtype='float')
        dictOut[names[n]] = {}
        dictOut[names[n]]['h'] = A[:,0]
        dictOut[names[n]]['k'] = A[:,1]
        dictOut[names[n]]['l'] = A[:,2]
        dictOut[names[n]]['FRobs'] = A[:,3]
        dictOut[names[n]]['FRcalc'] = A[:,4]
        dictOut[names[n]]['FRdiff'] = A[:,5]
        dictOut[names[n]]['Fnucl'] = A[:,6]
        dictOut[names[n]]['Fmag'] = A[:,7]
        dictOut[names[n]]['modFc'] = A[:,8]
        dictOut[names[n]]['scale'] = A[:,9]
        dictOut[names[n]]['weight'] = A[:,10]
    
    if params is not None:
        for n in range(len(names)):
            q = 2*np.pi/c.d_from_hkl(dictOut[names[n]]['h'],
                                    dictOut[names[n]]['k'],
                                    dictOut[names[n]]['l'],
                                    params)
                                    
            sort_indices = q.argsort()
            q = q[sort_indices]
            
            for entry in list(dictOut[names[n]].keys()):
                dictOut[names[n]][entry] = dictOut[names[n]][entry][sort_indices]
            
            dictOut[names[n]]['q'] = q
            dictOut[names[n]]['s'] = np.sqrt(2/dictOut[names[n]]['weight'])
            dictOut[names[n]]['fileName'] = names[n]
            dictOut[names[n]]['nVar'] = variables
    
    dictOut['variables'] = variables
    
    return dictOut

def calcR1(obs, calc, doPrint=False):
    
    N, D = 0, 0
    for n in range(len(obs)):
        N += abs(obs[n]-calc[n])
        D += abs(obs[n])

    R1 = N/D
    
    if doPrint:
        print('R1 = {}%'.format(round(R1*100,2)))
    return R1

def calcR2(obs, calc, doPrint=False):

    N, D = 0, 0
    for n in range(len(obs)):
        N += abs(obs[n]**2 - calc[n]**2)
        D += abs(obs[n]**2)
    
    R2 = N/D
    
    if doPrint:
        print('R2 = {}%'.format(round(R2*100,2)))
    return R2

def calcR3(obs, calc, w, doPrint=False):
    
    N, D = 0, 0
    for n in range(len(obs)):
        N += w[n]*abs(obs[n]-calc[n])
        D += w[n]*abs(obs[n])
        
    R3 = N/D
    
    if doPrint:
        print('R3 = {}%'.format(round(R3*100,2)))
    return R3
    
def calcR4(obs, calc, w, doPrint=False):

    N, D = 0, 0
    for n in range(len(obs)):
        N += w[n]*abs(obs[n]**2-calc[n]**2)
        D += w[n]*abs(obs[n]**2)

    R4 = N/D
    
    if doPrint:
        print('R4 = {}%'.format(round(R4*100,2)))
    return R4

def calcXsqrd(obs, calc, w, nParams, doPrint=False):

    N, D = 0, 0
    for n in range(len(obs)):
        N += w[n]*(obs[n]-calc[n])**2
    D += len(obs) - nParams

    Xsqrd = N/D
    
    if doPrint:
        print('X^2 = {:.3g}'.format(Xsqrd))
    return Xsqrd

def calcSWDsqrd(obs, calc, w, doPrint=False):

    SWDsqrd = 0
    for n in range(len(obs)):
        SWDsqrd += w[n]*(obs[n]-calc[n])**2
    
    if doPrint:
        print('SWD^2 = {:4.3e}'.format(SWDsqrd))
    return SWDsqrd