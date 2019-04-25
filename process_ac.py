import numpy as np
import os
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scientificConstants as sc
from scipy.optimize import curve_fit

globalPointMarker = 'o'
globalMarkerSize = 5
globalTextSize = 20

def calcTcolor(T, Tmin, Tmax):
    """
    Calculates the color that the corresponding line should be plotted with based on the
    idea that the coldest temperature should be completely blue RGB=(0,0,255) and the warmest
    temperature should be completely red RGB=(255,0,0)
    
    Input
    T: temperature of the curve that is being plotted
    D: dictionary containing all temperatures

    Output
    RGB: a tuple containing (R,G,B)-values for the color to be used
    """
    
    T18 = Tmin + 1/8*(Tmax-Tmin)
    T38 = Tmin + 3/8*(Tmax-Tmin)
    T58 = Tmin + 5/8*(Tmax-Tmin)
    T78 = Tmin + 7/8*(Tmax-Tmin)
    
    R, B, G = 0, 0, 0
    if T >= T78:
        p = (T-T78)/(Tmax-T78)
        R = 1-0.5*p
    elif T >= T58:
        p = (T-T58)/(T78-T58)
        R = 1
        G = 1-p
    elif T >= T38:
        p = (T-T38)/(T58-T38)
        R = p
        G = 1
        B = 1-p
    elif T >= T18:
        p = (T-T18)/(T38-T18)
        G = p
        B = 1
    else:
        p = (T-Tmin)/(T18-Tmin)
        B = 0.5+0.5*p
    
    return (R,G,B)
    
def numOfFreqRead(Temps):
    """
    Calculates the number of frequencies measured for each temperature in a PPMS data file
    
    Input
    Temps: a temperature column from a PPMS data file
    
    Output
    numOfFreq: an integer number of frequencies that has been measured
    """

    Temps = np.round(Temps, 1)*10
    Temps = np.array([int(x) for x in list(Temps)])
    bins = np.bincount(Temps)
    numOfFreq = int(np.max(bins))
    
    return numOfFreq

def readACDATAtoARRAY(fileName, dataOrigin='PPMS'):
    """Reads the output file from a PPMS ACMS measurement at AU and gives the result as a dictionary
    Input
    fileName: the name of the file to be read from
    numberofFreq: the number of frequency points for each temperature
    
    Output
    resultDict: a Python dictionary containing the specified below
    """
    
    colsToUse, rowsToSkip, delimiterToUse = 0, 0, 0
    if dataOrigin == 'PPMS':
        colsToUse = (2,3,4,5,8,9)
        rowsToSkip = 21
        delimiterToUse = ','
    
    f = open(fileName, 'r')
    names = f.readlines()
    f.close()
    
    names = [names[rowsToSkip-1].split(delimiterToUse)[n] for n in colsToUse]
    
    # Load the measured data into an array. Loading columns 'Temperature (K)', 'Magnetic Field (Oe)', 'Frequency (Hz)', 'Amplitude (Oe)', 'M' (emu)', 'M'' (emu)'
    D = np.loadtxt(fileName, delimiter=delimiterToUse, skiprows=rowsToSkip, usecols=colsToUse)
    
    # Reading the number of frequencies that have been measured
    numberofFreq = numOfFreqRead(D[:,0])
    #print('Splitting data according to {} frequencies per temperature'.format(numberofFreq))
    
    # Calculating how many different temperatures have been measured
    splitIntoArrays = int(len(D[:,0])/numberofFreq)
    
    # Splitting the original array into a list of arrays, each with its own temperature
    D = np.split(D, splitIntoArrays, axis=0)

    # Making a dictionary containing all array data. Each entry in the dictionary is a dictionary in itself,
    # which contains the data with the same temperature.
    resultDict = {}
    for d in D:
        T = np.mean(d[:,0]).round(1)
        temperatureDict = {}
        temperatureDict['T (K)'] = T
        for n in range(len(names)):
            temperatureDict[names[n]]=d[:,n]
        resultDict['{}K'.format(T)] = temperatureDict
        
    return resultDict

def calculateX(dict, molWeight, sampleMass, filmMass, setup='capsule'):
    """Adds lists of X' and X'' to the dictionary read from the AC data
    Input
    dict: the dictionary output from function readACDATAtoARRAY
    molWeight: molar weight of the sample in g/mol
    sampleMass: mass of the sample used in g
    filmMass: mass of the film used to wrap the sample in mg
    
    Output
    dict: the modified dictionary
    """
    
    Xd_capsule, Xd_film = 0, 0
    
    if setup == 'capsule':
        Xd_capsule = -1.8*10**-8 # unit: emu/Oe
        Xd_film = 6.47*10**-10 # unit: emu/(Oe*mg)
        
    else:
        print('Setup "{}" not recognized. Return.'.format(setup))
        return None
    
    Xd_sample = -0.6*10**-6
    
    # Modifying the individual dictionaries one at a time
    for key in dict.keys():
        # Read the content of the current dictionary into variable d
        d = dict[key]
        H = d['Amplitude (Oe)']
        H0 = d['Magnetic Field (Oe)']
        Mp = d["M' (emu)"]
        Mpp = d["M'' (emu)"]
        
        # Calculate X' and X'' from M' and M''
        Xp = (Mp - Xd_capsule*H - Xd_film*filmMass*H)*molWeight/(sampleMass*H) - Xd_sample*molWeight
        Xpp = Mpp/(sampleMass*H)*molWeight
        
        d["X' (emu/(Oe*mol))"] = np.array(Xp)
        d["X'' (emu/(Oe*mol))"] = np.array(Xpp)
    
    return dict    
    
def CCFITStartingGuesses(D):
    """
    Reads the dictionary prepared by readACDATAtoARRAY and estimates starting parameters for CCFit based on the CCFit manual
    
    Input
    D: dictionary containing PPMS data in the form as returned by readACDATAtoARRAY
    
    Output:
    res: tuple containing (Xs, Xt, tau, alpha)
    """
    
    # Reading the lowest measured temperature and getting the result as a string
    T = sorted([float(x[0:-1]) for x in list(D.keys())])[0] # the lowest temperature
    Ts = '{}K'.format(T) # the lowest temperature as a string
    
    # Extracting the data from the lowest temperature measurement
    v = D[Ts]['Frequency (Hz)']
    Xp = D[Ts]["X' (emu/(Oe*mol))"]
    Xpp = D[Ts]["X'' (emu/(Oe*mol))"]
    iXppMax = np.argmax(np.array(Xpp))
    tau = 1/(2*np.pi*v[iXppMax])
    Xs = 0
    Xt = Xp[0]
    alpha = 0.1
    
    res = (Xs, Xt, tau, alpha)

    return res

def prepareCCFITinput(D, processes=1):
    """Writes an input file for CCfit from the Chilton group"""
    
     # Sorting the temperatures for printing of the cc-fit input file
    Ts = list(D.keys())
    Ts = [float(t[0:-1]) for t in Ts]
    Ts.sort()
    
    # Opening input file
    f = open('ccinput.dat', 'w')
    
    # Reading the number of temperatures and frequencies from the data and writes them to file head
    temperatures = len(Ts)
    frequencies = len(D[list(D.keys())[0]]['Frequency (Hz)'])
    f.write('{} {} {}\n'.format(processes, temperatures, frequencies))
    
    # Inputting starting parameters into file head
    f.write('{} {} {} {}\n'.format(*CCFITStartingGuesses(D)))
    
    # Writing the frequency and susceptibility data to the file according to the manual specifications
    for n in range(len(D['{}K'.format(Ts[0])]['Frequency (Hz)'])):
        s = '{} '.format(D['{}K'.format(Ts[0])]['Frequency (Hz)'][n])
        for T in Ts:
            Xp = D['{}K'.format(T)]["X' (emu/(Oe*mol))"][n]
            Xpp = D['{}K'.format(T)]["X'' (emu/(Oe*mol))"][n]
            s += '{} {} '.format(Xp, Xpp)
        s += '\n'
        f.write(s)
    
    f.close()
    
def runCCFIT(D):
    """
    Runs the program CCFit in the folder of the script
    """
    
    if 'ccinput.dat' not in os.listdir():
        prepareCCFITinput(D)
        print('Input file for cc-fit not found. Preparing file using 1 relaxation process.')
    
    CCFIT = Popen('cc-fit ccinput.dat', stdin=PIPE, stdout=PIPE, stderr=PIPE)
    CCFIT.stdin.write(b'\n')
    out, err = CCFIT.communicate()
    
    sortedTemps = ['{}K'.format(T) for T in sorted([float(x[0:-1]) for x in list(D.keys())])]
    data = readCCout()
    
    n = 0
    while n < len(sortedTemps):
        D[sortedTemps[n]]['ccXs'] = data[n,0]
        D[sortedTemps[n]]['ccXt'] = data[n,1]
        D[sortedTemps[n]]['cctau'] = data[n,2]
        D[sortedTemps[n]]['ccalpha'] = data[n,3]
        D[sortedTemps[n]]['ccresidual'] = data[n,4]
        n += 1
    
    return D
    
def readCCout():
    """
    Reads the .out-file containing the fitted parameters from CCFit
    
    Output
    data: a n-by-5 array containing the fitted parameters and residuals from the .out file
    """

    inputName = 'ccinput.dat'
    ccextension = '_cc-fit.out'
    
    data = np.loadtxt(inputName+ccextension, skiprows = 13)
    
    return data
    
def Xp_(v, Xs, Xt, tau, alpha):
    """
    Calculates the function X' [chi prime] as specified in Molecular Nanomagnets eq. 3.27
    
    Input:
    v: frequency of the AC field
    Xs: adiabatic limit of susceptibility
    Xt: isothermal limit of susceptibility
    tau: relaxation time of the system
    alpha: width of relaxation time distribution
    
    Output
    Xp: the function value at the specified frequency
    """

    w = 2*np.pi*v
    
    Xp = Xs + (Xt - Xs)*(1 + (w*tau)**(1-alpha)*np.sin(np.pi*alpha/2))/(1 + 2*(w*tau)**(1-alpha)*np.sin(np.pi*alpha/2) + (w*tau)**(2-2*alpha))
    
    return Xp

def Xpp_(v, Xs, Xt, tau, alpha):
    """
    Calculates the function X'' [chi double-prime] as specified in Molecular Nanomagnets eq. 3.27
    
    Input:
    v: frequency of the AC field
    Xs: adiabatic limit of susceptibility
    Xt: isothermal limit of susceptibility
    tau: relaxation time of the system
    alpha: width of relaxation time distribution
    
    Output
    Xpp: the function value at the specified frequency
    """

    w = 2*np.pi*v
    
    Xpp = (Xt - Xs)*(w*tau)**(1-alpha)*np.cos(np.pi*alpha/2)/(1 + 2*(w*tau)**(1-alpha)*np.sin(np.pi*alpha/2) + (w*tau)**(2-2*alpha))

    return Xpp

def plotXppvsFreq(D, addFit=False):
    """Plots frequency vs. X'' for each temperature"""
    
    Ts = sorted([float(x[0:-1]) for x in list(D.keys())])
    Tmax = Ts[-1]
    Tmin = Ts[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for key in D.keys():
        d = D[key]
        if addFit:
            d = D[key]
            Xs = d['ccXs']
            Xt = d['ccXt']
            tau = d['cctau']
            alpha = d['ccalpha']
            v = np.linspace(d['Frequency (Hz)'][0], d['Frequency (Hz)'][-1], 1000)
            ax.semilogx(v, Xpp_(v, Xs, Xt, tau, alpha), c=calcTcolor(d['T (K)'], Tmin, Tmax))
            ax.semilogx(d['Frequency (Hz)'], d["X'' (emu/(Oe*mol))"], marker=globalPointMarker,
                                                                      markerfacecolor='None',
                                                                      markersize=globalMarkerSize,
                                                                      linestyle='None',
                                                                      markeredgecolor='k',
                                                                      linewidth=1)
        
        else:
            ax.semilogx(d['Frequency (Hz)'], d["X'' (emu/(Oe*mol))"], c=calcTcolor(d['T (K)'], Tmin, Tmax), marker='o', linestyle='None')
    
    large_T = mpatches.Patch(color=(0.5,0,0), label='{}K'.format(Tmax))
    small_T = mpatches.Patch(color=(0,0,0.5), label='{}K'.format(Tmin))

    ax.tick_params(labelsize=globalTextSize)    
    ax.set_xlabel(r"$\nu$ [$Hz$]", fontsize=globalTextSize)
    ax.set_ylabel(r"$\chi$'' [$\frac{emu}{Oe  mol}$]", fontsize=globalTextSize)
    ax.legend(handles=[small_T, large_T], fontsize=globalTextSize)
    
    return fig
    
def ColeColePlot(D, addFit=False):
    """Plots a Cole-Cole plot of the measured data"""

    Ts = sorted([float(x[0:-1]) for x in list(D.keys())])
    Tmax = Ts[-1]
    Tmin = Ts[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for key in D.keys():
        d = D[key]
        if addFit:
            d = D[key]
            Xs = d['ccXs']
            Xt = d['ccXt']
            tau = d['cctau']
            alpha = d['ccalpha']
            v = np.linspace(d['Frequency (Hz)'][0], d['Frequency (Hz)'][-1], 1000)
            ax.plot(Xp_(v, Xs, Xt, tau, alpha), Xpp_(v, Xs, Xt, tau, alpha), c=calcTcolor(d['T (K)'], Tmin, Tmax), linestyle='-')
            ax.plot(d["X' (emu/(Oe*mol))"], d["X'' (emu/(Oe*mol))"], marker=globalPointMarker,
                                                                     markerfacecolor='None',
                                                                     markersize=globalMarkerSize,
                                                                     linestyle='None',
                                                                     markeredgecolor='k',
                                                                     linewidth=1)
    
        else:
            ax.plot(d["X' (emu/(Oe*mol))"], d["X'' (emu/(Oe*mol))"], c=calcTcolor(d['T (K)'], Tmin, Tmax), marker='o', linestyle='None')

    large_T = mpatches.Patch(color=(0.5,0,0), label='{}K'.format(Tmax))
    small_T = mpatches.Patch(color=(0,0,0.5), label='{}K'.format(Tmin))
    
    ax.tick_params(labelsize=globalTextSize)
    ax.set_xlabel(r"$\chi$' [$\frac{emu}{Oe  mol}$]", fontsize=globalTextSize)
    ax.set_ylabel(r"$\chi$'' [$\frac{emu}{Oe  mol}$]", fontsize=globalTextSize)
    ax.legend(handles=[small_T, large_T], fontsize=globalTextSize)
        
    return fig

def tauvsTempPlot(D):
    """
    Plots the relaxation time vs. temperature for the measured data.
    
    Input
    D: the dictionary that has been modified by CCfit to incorporate tau
    
    Output:
    res: a tuple containing the arrays (T, tau)
    """
    
    T = sorted([float(x[0:-1]) for x in list(D.keys())])
    tau = []
    
    for t in T:
        tau.append(D['{}K'.format(t)]['cctau'])
    
    T, tau = np.array(T), np.array(tau)
    
    plt.plot(1/T, np.log(tau), 'bo')
    plt.show()
    
    return tuple([T, tau])

def readXppvsT(D):

    Ts = []
    for key in D.keys():
        Ts.append(float(key[:-1]))
    Ts = sorted(Ts)
    
    Q = {}
    W = D['{}K'.format(Ts[0])]['Frequency (Hz)']
    for w in W:
        Q['{}Hz'.format(round(w,2))] = {'Temp (K)':[], 'Xp':[], 'Xpp':[]}
    for T in Ts:
        w = D['{}K'.format(T)]['Frequency (Hz)']
        Xp = D['{}K'.format(T)]["X' (emu/(Oe*mol))"]
        Xpp = D['{}K'.format(T)]["X'' (emu/(Oe*mol))"]
        for n in range(len(w)):
            Q['{}Hz'.format(round(w[n],2))]['Temp (K)'].append(T)
            Q['{}Hz'.format(round(w[n],2))]['Xp'].append(Xp[n])
            Q['{}Hz'.format(round(w[n],2))]['Xpp'].append(Xpp[n])
    
    return Q

def makeOrbachfit(fileName, minT, maxT, U_guess, t0_guess, U_unit='K'):
    """
    Reads the (T,tau)-data from <fileName> and does a fitting to an
    Orbach-process.
    
    Input
    fileName: name of the file in the format printed by makeTvsTauFile
    minT: minimum temperature used for the fit
    maxT: maximum temperature used for the fit
    t0_guess: guess for the t0-parameter in the Orbach model
    U_guess: guess for the U-parameter in the Orbach-model
    U_guess_unit: unit of the guess given for U-parameter (default: 'K' (Kelvin))
    
    Output
    P: parameters and their uncertainties from the modeling
    """
    
    D = np.loadtxt(fileName, skiprows=1)
    T = D[:,0]
    tau = D[:,1]
    
    # Calculating the values to use for the fit
    T_used = np.array([x for x in T if x>=minT and x<=maxT])
    T_not_used = np.array([x for x in T if x<minT or x>maxT])
    
    # Calculating the values not to use for the fit
    tau_used = np.array([tau[np.where(T==x)[0][0]] for x in T_used])
    tau_not_used = np.array([tau[np.where(T==x)[0][0]] for x in T_not_used])
    
    # Recalculating U_guess
    if U_unit=='K':
        U_guess = U_guess*sc.kB
    elif U_unit=='J':
        pass
    else:
        print('Wrong unit given. Ending.')
        return
    
    # Fitting
    A = curve_fit(OrbachTau, T_used, tau_used, [U_guess, t0_guess])
    
    U_fitted = A[0][0]
    t0_fitted = A[0][1]
    
    uncertainties = np.sqrt(np.diag(A[1]))
    U_sigma = uncertainties[0]
    t0_sigma = uncertainties[1]
    
    # Recalculating U_fitted
    if U_unit=='K':
        U_fitted = U_fitted/sc.kB
        U_sigma = U_sigma/sc.kB
    elif U_unit=='J':
        pass
    
    plt.title('no_linear')
    plt.plot(1/T_used, np.log(OrbachTau(T_used, A[0][0], A[0][1])))
    plt.plot(1/T_used, np.log(tau_used),'k*')
    plt.plot(1/T_not_used, np.log(tau_not_used), 'r*')
    
    P = (U_fitted, U_sigma, t0_fitted, t0_sigma)
    print('no_linear')
    print('Ueff = {:3.2f} +/- {:3.2f}\nt0 = {:3.2f} +/- {:3.2f}'.format(P[0], P[1], P[2], P[3]))
    plt.show()
    
    return P

def OrbachTau(T, Ueff, t0):
    """
    Returns the Orbach relaxation time at a given temperature according to the parameters Ueff and t0
    
    Input
    T: temperature list or float
    Ueff: effective energy barrier for Orbach relaxation in joules
    t0: the maximal Orbach relaxation rate
    
    Output
    tau: a list or float giving the Orbach relaxation times
    """
    
    tau = t0*np.exp(Ueff/(sc.kB*T))
    
    return tau

def RamanTau(T, Ct, n):
    """
    Returns the Raman relaxation time at a given temperature according to the parameters Ct and n
    
    Input
    T: temperature list or float
    Ct: characteristic constant for Raman relaxation
    n: the Raman exponent of the relaxation
    
    Output
    tau: a list or float giving the Raman relaxation times
    """

    tau = Ct*T**(-n)
    
    return tau
    
def QTTau(T, tQT):
    """
    Returns the quantum tunneling relaxation time at a given temperature according to the parameter tQT
    
    Input
    T: temperature list or float
    tQT: quantum tunneling relaxation time
    """
    
    tau = tQT
    
    return tau
    
def getParameterGuesses(T, tau):
    """
    Calculates guesses for optimal fitting parameters to begin the fit
    
    Input
    T: temperature array
    tau: relaxation time array
    
    Output
    guess: dictionary of guessed parameters for the relaxation functions
    """
    
    # Obtaining points for Orbach parameter guessing
    T1, tau1 = T[-1], tau[-1]
    T2, tau2 = T[-2], tau[-2]
    
    # Calculating guesses for Orbach parameters
    Ueff_guess = (np.log(tau2) - np.log(tau1))/(1/T2-1/T1)*sc.kB
    
    t0_guess = tau1*np.exp(-Ueff_guess/(sc.kB*T1))
    
    # Obtaining points for Raman parameter guessing
    l = len(T)
    index1, index2 = 0,0
    if l/2 % 1 == 0:
        index1, index2 = int(l/2), int(l/2+1)
    else:
        index1, index2 = int(l/2-0.5), int(l/2+0.5)
    
    T1, tau1 = T[index1], tau[index1]
    T2, tau2 = T[index2], tau[index2]    
    
    # Calculating guesses for Raman parameters
    n_guess = (np.log(tau1) - np.log(tau2))/(np.log(T2) - np.log(T1))
    
    Cr_guess = tau1*T1**n_guess
    
    # Extracting guess for QT parameter
    tQT_guess = tau[0]
    
    guess = {'Ueff': Ueff_guess,
            't0': t0_guess,
            'n': n_guess,
            'Cr': Cr_guess,
            'tQT': tQT_guess}
    
    return guess
    
def makeTvsTauFile(D, fileName, saveFile=True):
    """
    Prints the data from the dictionary D into a file containing 2 columns of
    temperature data and relaxation time data
    
    Input
    D: dictionary as returned from runCCFIT
    fileName: name to give to the new file that the function creates
    
    Output
    T: sorted array of measured temperatures
    tau: array matching the sorting of T with fitted relaxation times
    """
    
    T = sorted([float(x[:-1]) for x in D.keys()])
    tau = [D['{}K'.format(x)]['cctau'] for x in T]
    
    if saveFile:
        f = open(fileName, 'w')
        f.write('{:>6s}{:>15s}\n'.format('T', 'tau'))
    
        for n in range(len(T)):
            f.write('{:6.3f}{:15.6e}\n'.format(T[n], tau[n]))
        
        f.close()
        
        print('Data written to file {}'.format(fileName))
    
    return np.array(T), np.array(tau)

def _QT(T, tQT):
    """
    Basic function for calculating relaxation time due to
    quantum tunneling
    
    Input
    T: temperature for the calculation
    tQT: characteristic time for quantum tunneling
    
    Output
    tau: relaxation time due to quantum tunneling
    """
    
    tau = tQT
    
    return tau

def _R(T, Cr, n):
    """
    Basic function for calculating relaxation time due to
    the Raman mechanism
    
    Input
    T: temperature for the calculation
    Cr: Raman pre-factor
    n: Raman exponent
    
    Output
    tau: relaxation time due to the Raman mechanism
    """
    
    tau = Cr*T**-n

    return tau
    
def _O(T, t0, Ueff):
    """
    Basic function for calculating relaxation time due to
    the Orbach relaxation mechanism
    
    Input
    T: temperature for the calculation
    t0: characteristic time for quantum tunneling
    Ueff: effective energy barrier to thermal relaxation
    
    Output
    tau: relaxation time due to the Orbach mechanism
    """
    
    tau = t0*np.exp(Ueff/(sc.kB*T))
    
    return tau
    
def _QT_log(T, tQT):
    """
    Wrapper function to function _QT that computes the logarithmic
    relaxation time due to quantum tunneling.
    
    See help(_QT) for more
    """
    
    return np.log(_QT(T, tQT))

def _R_log(T, Cr, n):
    """
    Wrapper function to function _R that computes the logarithmic
    relaxation time due to the Raman mechanism.
    
    See help(_R) for more
    """    
    
    return np.log(_R(T, Cr, n))
    
def _O_log(T, t0, Ueff):
    """
    Wrapper function to function _O that computes the logarithmic
    relaxation time due to the Orbach mechanism.
    
    See help(_O) for more
    """
    
    return np.log(_O(T, t0, Ueff))
    
def _QTR(T, tQT, Cr, n):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism and the Raman mechanism
    
    See help(_QT) and help(_R) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_R(T, Cr, n)
    
    tau = 1/w
    
    return np.log(tau)

def _QTO(T, tQT, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism and the Orbach mechanism
    
    See help(_QT) and help(_O) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)
    
def _RO(T, Cr, n, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a Raman
    mechanism and the Orbach mechanism
    
    See help(_R) and help(_O) for more
    """
    
    w = 1/_R(T, Cr, n) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)

def _QTRO(T, tQT, Cr, n, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism, the Raman mechanism and the Orbach mechanism
    
    See help(_QT), help(_R) and help(_O) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_R(T, Cr, n) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)

def getStartParams(guess, fitType='QTRO'):
    
    p0 = 0
    if fitType=='QT':
        p0 = [guess['tQT']]
    elif fitType=='R':
        p0 = [guess['Cr'], guess['n']]
    elif fitType=='O':
        p0 = [guess['t0'], guess['Ueff']]
    elif fitType=='QTR':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='R')
    elif fitType=='QTO':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='O')
    elif fitType=='RO':
        p0 = getStartParams(guess, fitType='R') + getStartParams(guess, fitType='O')
    elif fitType=='QTRO':
        p0 = [guess['tQT'], guess['Cr'], guess['n'], guess['t0'], guess['Ueff']]
    else:
        print('fitType parameter did not correspond to any correct one')
        
    return p0

def getFittingFunction(fitType='QTRO'):
    
    f = 0
    if fitType=='QT':
        f = _QT_log
    elif fitType=='R':
        f = _R_log
    elif fitType=='O':
        f = _O_log
    elif fitType=='QTR':
        f = _QTR
    elif fitType=='QTO':
        f = _QTO
    elif fitType=='RO':
        f = _RO
    elif fitType=='QTRO':
        f = _QTRO
    else:
        print('fitType parameter did not correspond to any correct one')
        
    return f
    
def readPopt(popt, pcov, fitType='QTRO'):

    p_fit = 0
    if fitType=='QT':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT']}
    elif fitType=='R':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n']}
    elif fitType=='O':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['t0', 'Ueff']}
    elif fitType=='QTR':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n']}
    elif fitType=='QTO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 't0', 'Ueff']}
    elif fitType=='RO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n', 't0', 'Ueff']}
    elif fitType=='QTRO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n', 't0', 'Ueff']}
    
    return p_fit
    
def fitRelaxation(fileName, tempRange, fitType='QTRO'):
    """
    Documentation to come
    """
    
    D = np.loadtxt(fileName, skiprows=1)
    T = D[:,0]
    tau = D[:,1]
    minT = tempRange[0]
    maxT = tempRange[1]
    
    # Calculating the values to use for the fit
    T_used = np.array([x for x in T if x>=minT and x<=maxT])
    T_not_used = np.array([x for x in T if x<minT or x>maxT])
    
    T_used_space = np.linspace(T_used[0], T_used[-1], 1000)
    
    # Calculating the values not to use for the fit
    tau_used = np.array([tau[np.where(T==x)[0][0]] for x in T_used])
    tau_not_used = np.array([tau[np.where(T==x)[0][0]] for x in T_not_used])
    
    # Obtaining automated guesses for fitting parameters
    guess = getParameterGuesses(T, tau)
    f = getFittingFunction(fitType=fitType)
    p0 = getStartParams(guess, fitType=fitType)
    
    # Making a fit to the selected data
    popt, pcov = curve_fit(f, T_used, np.log(tau_used), p0)
    p_fit = readPopt(popt, pcov, fitType=fitType)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(1/T_used_space, np.ones(T_used_space.shape)*f(T_used_space, *p_fit['params']), 'b-', label='Fitted')
    ax.plot(1/T_used, np.log(tau_used), marker=globalPointMarker,
                                        markersize=globalMarkerSize,
                                        markeredgecolor='k',
                                        markerfacecolor='None',
                                        linestyle='None')
                                        
    ax.plot(1/T_not_used, np.log(tau_not_used), marker=globalPointMarker,
                                                markersize=globalMarkerSize,
                                                markeredgecolor='r',
                                                markerfacecolor='None',
                                                linestyle='None')
    
    ax.set_xlabel(r'$\frac{1}{T}$ [$K^{-1}$]', fontsize=globalTextSize)
    ax.set_ylabel(r'$\ln{\tau}$ [$\ln{s}$]', fontsize=globalTextSize)
    ax.tick_params(labelsize=globalTextSize)
    ax.legend(fontsize=globalTextSize)
    
    return fig, p_fit
    
def printFittedParams(p_fit):
    """
    Prints the parameters that are in the dictionary p_fit
    with the names and uncertainties
    """
    
    print('\nThe fitted parameters are:')
    for n in range(len(p_fit['params'])):
        print('{} = {} +/- {}'.format(p_fit['quantities'][n], p_fit['params'][n], p_fit['sigmas'][n]))

def readPFITinOrder(p_fit, plotType='O'):

    p = []
    if plotType=='QT':
        tQT = p_fit['params'][p_fit['quantities'].index('tQT')]
        p = [tQT]
    elif plotType=='R':
        Cr = p_fit['params'][p_fit['quantities'].index('Cr')]
        n = p_fit['params'][p_fit['quantities'].index('n')]
        p = [Cr, n]
    elif plotType=='O':
        t0 = p_fit['params'][p_fit['quantities'].index('t0')]
        Ueff = p_fit['params'][p_fit['quantities'].index('Ueff')]
        p = [t0, Ueff]
    elif plotType=='QTR':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R')
    elif plotType=='QTO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='RO':
        p = readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='QTRO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    
    return p

def readTvsTauFile(fileName):
    
    D = np.loadtxt(fileName, skiprows=1)
    
    T = D[:,0]
    tau = D[:,1]
    
    return T, tau
    
def addPartialModel(fig, T, tau, Tmin, Tmax, p_fit, plotType='O'):
    
    ax = fig.get_axes()[0]
    
    guess = getParameterGuesses(T, tau)
    f = getFittingFunction(fitType=plotType)
    p = readPFITinOrder(p_fit, plotType=plotType)
    
    T_space = np.linspace(Tmin, Tmax, 1000)
    ax.plot(1/T_space, np.ones(T_space.shape)*f(T_space, *p), label=plotType)
    ax.legend(fontsize=globalTextSize)

def readTwoColumnCSV(fileName, newFile, saveFile=False):

    D = np.loadtxt(fileName, delimiter=';')
    D = np.split(D, len(D[:,0])/4)
    
    T_recip, lntau = [], []
    for ary in D:
        ary = np.mean(ary, axis=0)
        T_recip.append(ary[0])
        lntau.append(ary[1])
    
    T = 1/np.array(T_recip)
    tau = np.exp(np.array(lntau))
    
    sort_indices = T.argsort()
    T = T[sort_indices]
    tau = tau[sort_indices]
    
    if saveFile:
        f = open('{}.dat'.format(fileName.split('.')[1]), 'w')
        for n in range(len(T)):
            f.write('{} {}\n'.format(T[n], tau[n]))
        f.close()
    
    return T, tau

def printUeffInKelvin(p_fit):
    """
    Reads the dictionary p_fit to see, whether an Orbach fit was made.
    If the fit was made, the effective barrier is printed in Kelvin.
    If the fit was not made, a message is printed telling that.
    """
    
    quantities = p_fit['quantities']
    
    if 'Ueff' in quantities:
        index = quantities.index('Ueff')
        Ueffval = p_fit['params'][index]/sc.kB
        Ueffsigma = p_fit['sigmas'][index]/sc.kB
        
        print('Ueff = {}K +/- {}K'.format(Ueffval, Ueffsigma))
    
    else:
        print('A Ueff-parameter has not been produced by the fit.')

def plotSusceptibility(D, T_plot, type='Xpp'):

    T = np.array(sorted([float(x[:-1]) for x in D.keys()]))
    index = np.argmin(np.abs(T-T_plot))
    T_to_use = T[index]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x,y = 0,0
    if type == 'Xpp':
        y = D['{}K'.format(T_to_use)]["X'' (emu/(Oe*mol))"]
        type = "''"
    elif type == 'Xp':
        y = D['{}K'.format(T_to_use)]["X' (emu/(Oe*mol))"]
        type = "'"
    x = D['{}K'.format(T_to_use)]['Frequency (Hz)']
    
    ax.semilogx(x,y, marker=globalPointMarker,
                     linestyle='None',
                     markerfacecolor=calcTcolor(T_to_use, 4.0, 25),
                     markeredgecolor='None',
                     label='{}K'.format(T_to_use))
    ax.set_title(r'$\chi${} vs frequency at {}K'.format(type, T_to_use), fontsize=globalTextSize)
    ax.set_xlabel(r"$\nu$ [$Hz$]", fontsize=globalTextSize)
    ax.set_ylabel(r'$\chi$'+'{}'.format(type)+r'[$\frac{emu}{Oe  mol}$]', fontsize=globalTextSize)
    ax.legend()
    
    return fig
    
if __name__ == '__main__':
    
    """
    PUT YOUR VALUES HERE
    """
    # Name of the file printed by the PPMS
    fileName = 'data\\dy_fod\\20180430DyfodOPy-ac-0.dat'
    
    # Molar weight of the sample
    molWeight = 2286.18 # in g/mol
    
    # Mass of the sample that was used (in g)
    sampleMass = 32.8*10**-3 # in g
    
    # Mass of the Parafilm that was used (in mg)
    filmMass = 13.7

    # The temperature range in which the fit should be performed
    tempRange = [11.5, 20]

    # Whether to plot the result of CCFit on top of measured data
    plotFit = True

    # Whether to save the TvsTau file or not
    saveFile = True

    # Name of the TvsTau-file, if it is to be made
    TvsTaufileName = 'Dy_FOD_ac0Oe.dat'

    # Which type of fit should be made? (In order QT, R, O. For example to use quantum tunneling and Raman, put in 'QTR'. The opposite, namely 'RQT' will not work)
    fitToMake = 'RO'
    
    D = readACDATAtoARRAY(fileName, dataOrigin='PPMS')
    D = calculateX(D, molWeight, sampleMass, filmMass, setup='capsule')
    prepareCCFITinput(D)
    D = runCCFIT(D)
    T, tau = makeTvsTauFile(D, TvsTaufileName, saveFile=saveFile)
    guess = getParameterGuesses(T, tau)
    
    fig1 = ColeColePlot(D, addFit=plotFit)
    fig2 = plotXppvsFreq(D, addFit=plotFit)
    fig3, p_fit = fitRelaxation(TvsTaufileName, tempRange, fitType=fitToMake)
    
    # The function addPartialModel can be used to add a plot showing the T vs tau simulation for a single process.
    # Use it as many times as deemed appropriate.
    # Parameters (in order):
    # <figure to plot on>,
    # <a list of temperatures>, 
    # <a list of relaxation times>, 
    # <minimum temperature for simulation>, 
    # <maximum temperature for simulation>, 
    # <fitted parameters from the fitRelaxation-function>, 
    # <which type of plot to add (possibilities: 'QT' or 'R' or 'O'>
    
    # For example, to add a simulation of a Raman process to figure 3 between
    # 7 K and 18 K, use the line below
    # addPartialModel(fig3, T, tau, 7, 18, p_fit, plotType='R')
    
    # For example, to add a simulation of an Orbach process to figure 3 between
    # 15 K and 25 K, use the line below
    # addPartialModel(fig3, T, tau, 15, 25, p_fit, plotType='O')
    
    # KEEP THESE METHOD CALLS
    printFittedParams(p_fit)
    printUeffInKelvin(p_fit)
    ax = fig3.get_axes()[0]
    ax.legend()
    plt.show()











