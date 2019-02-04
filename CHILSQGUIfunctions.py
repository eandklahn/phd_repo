import os
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import calculateCHILSQRs as cc
import numpy as np
import time

def importCRYfile():

    cryFile = [f for f in os.listdir() if f.endswith('.cry')][0]
    
    f = open(cryFile, 'r')
    data = f.readlines()
    f.close()
    
    chiNumbers = []
    chiIons = []
    for n in range(len(data)):
        line = data[n]
        if line.split()[0] == 'Q' and line.split()[2] == 'CHI':
            chiIons.append(line.split()[1])
            chiNumbers.append(n)
    
    return [data, chiNumbers, cryFile, chiIons]
    

def changeASPs(D, atom, CH11, CH22, CH33, CH23, CH31, CH12):

    data = D[0]
    chiNumbers = D[1]
    cryFile = D[2]
    N = 0
    for n in chiNumbers:
        line = data[n].split()
        if line[1] == atom:
            N = n
    
    CH11 = float(CH11)
    CH22 = float(CH22)
    CH33 = float(CH33)
    CH23 = float(CH23)
    CH31 = float(CH31)
    CH12 = float(CH12)
    
    newData = data[:N]
    newData.append('Q {} CHI {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f}\n'.format(
                                                                                atom, CH11, CH22, CH33, CH23, CH31, CH12))
    newData.extend(data[N+1:])
    
    f = open(cryFile, 'w')
    for line in newData:
        f.write(line)
    f.close()
    
    return [newData, chiNumbers, cryFile, D[3]]    

def getASPs(data, chiNumbers, atom):
    
    for n in chiNumbers:
        line = data[n].split()
        if line[1] == atom:
            return (float(line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7]), float(line[8]))
    
def runCHILSQ(cryFile, runNumber):
    
    if 'chilsq.lis' in os.listdir():
        print('Delete chilsq.lis first')
        return
    
    extFiles = [f for f in os.listdir() if f.endswith('.ext')]
    
    for file in extFiles:
        string = '{}\n{} 5 5 2\n\n\n'.format(cryFile, file)
        
        chilsq = Popen('chilsq', stdin=PIPE, stdout=PIPE, stderr=PIPE)
        chilsq.stdin.write(bytes(string, 'utf8'))
        chilsq.communicate()

        os.rename('chilsq.lis', 'chilsq{}.lis'.format(runNumber))
        runNumber += 1
    
    return runNumber
    
def handleOldListingsOnStartup():
    
    lisFiles = [f for f in os.listdir(os.getcwd()) if f.endswith('.lis')]

    for f in lisFiles:
        os.rename(f, os.getcwd()+'\\oldListings\\'+f)
    
def getRunNumber():
    
    L = []
    for l in os.listdir(os.getcwd()+'\oldListings'):
        l = l.split('.')[0][6:]
        if l == '':
            continue
        else:
            l = int(l)
            L.append(l)

    if len(L) == 0:
        return 1
    else:
        L.sort()
        return int(L[-1]+1)
    
def createListingsFolder():

    if 'oldListings' in os.listdir():
        return getRunNumber()
    else:
        os.mkdir('oldListings')
        return 1

def calculateXsqrdwithCHILSQ(p0):
    """
    Function to calculate Xsqrd from a set of parameters.
    This function is obsolete, as the functionality has been put in a method within
    the CHILSQfit-class. In there, 'Dy1' is for example no longer hard coded
    """
    
    D = importCRYfile()
    
    # Remember that 'Dy1' is hardcoded at the moment - change that!
    D = changeASPs(D, 'Dy1', *p0)
    
    runCHILSQ(D[2], 0)
    
    lisFiles = [f for f in os.listdir() if f.endswith('.lis')]
    
    allFRobs = []
    allFRcalc = []
    allweights = []
    
    for l in lisFiles:
        D = cc.readCHILSQ_FRs(l)
        rflName = list(D.keys())[0]
        dict = D[rflName]
        FRobs, FRcalc, weight = (list(dict['FRobs']),
                                 list(dict['FRcalc']),
                                 list(dict['weight']))
        
        allFRobs += FRobs
        allFRcalc += FRcalc
        allweights += weight
        
        #print(cc.calcXsqrd(np.array(FRobs), np.array(FRcalc), np.array(weight), 0))
    
    for l in lisFiles:
        os.remove(l)
    
    allFRobs = np.array(allFRobs)
    allFRcalc = np.array(allFRcalc)
    allweights = np.array(allweights)
    
    X_sqrd = cc.calcXsqrd(allFRobs, allFRcalc, allweights, 0)
    
    return X_sqrd

def CHILSQcostFunction(p):
    """
    Written as a wrapper function for calculateXsqrdwithCHILSQ
    to calculate a cost function from Xsqrd.
    
    This function is obsolete. We can minimize Xsqrd directly
    without the need for another cost function.
    """
    
    X_sqrd = calculateXsqrdwithCHILSQ(p)
    
    print('Current value of X squared: {:6.3f} w params {:4.3f}, {:4.3f}, {:4.3f}, {:4.3f}, {:4.3f}, {:4.3f}'.format(X_sqrd, *p))
    
    return X_sqrd
    
if __name__ == '__main__':

    t_start = time.time()

    lisFiles = [f for f in os.listdir() if f.endswith('.lis')]
    for l in lisFiles:
        os.remove(l)
    
    #p0 = [0.1,0.1,0.1,0.1,0.1,0.1]
    
    p_fit = minimize(CHILSQcostFunction, p0, method='Powell')
    
    t_end = time.time()
    
    print('Finished in {:6.2f} seconds'.format(t_end-t_start))
    
    p = p_fit.x
    
    CH11 = p[0]
    CH22 = p[1]
    CH33 = p[2]
    CH23 = p[3]
    CH31 = p[4]
    CH12 = p[5]
    
    X = np.array([[CH11, CH12, CH31],
                  [CH12, CH22, CH23],
                  [CH31, CH23, CH33]])
    
    T, V = np.linalg.eig(X)
    
    print("""Obtained chi-tensor
{}
with eigenvalues
{}""".format(X, T))
    