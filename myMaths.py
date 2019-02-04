import matplotlib.pyplot as plt
import numpy as np

def sqrt(x):
    """Finds the square root of the number x by using Newtons method"""
    
    guess = 1
    
    while abs(x/guess - guess) > 10**-12:
        guess = (guess + x/guess)/2
        
    return guess

def isPrime(N,Ps):
    """Checks if a number N is prime given a list of primes Ps where the primes are all smaller than sqrt(N)"""
    
    n = 0
    while Ps[n] <= sqrt(N):
        if N % Ps[n] == 0:
            return False
        else:
            n += 1
    return True

def prime_pi(x):
    """Calculates the number of primes less than or equal to x"""
    
    Ps = [2,3]
    if x < 4:
        return 2
        if x < 3:
            return 1
            if x < 2:
                return 0
    else:
        for N in range(4,x+1):
            if isPrime(N,Ps):
                Ps.append(N)
    
    return len(Ps)

def gaussian(x,A,x0,s):
    """Calculates a Gaussian function
    x: plotting range
    A: amplitude
    x0: x-offset
    s: uncertainty
    """
    
    return A*np.exp(-(x-x0)**2/(2*s**2))

def lorentzian(x,A,x0,s):
    """Calculates a Lorentzian function
    x: plotting range
    A: amplitude
    x0: x-offset
    s: uncertainty
    """
    
    return A*(s**2/((x-x0)**2+s**2))
    
def pseudo_voigt(x,mu,x0,I,sg,sl):
    """Calculates a Pseudo-Voigt function
    x: plotting range
    A: amplitude
    x0: x-offset
    s: uncertainty
    """
    assert mu > 0 and mu < 1
    
    return mu*lorentzian(x,I,x0,sl)+(1-mu)*gaussian(x,I,x0,sg)
    
def division(a, b, e=1e-15):
    """
    Divides a by b: a/b
    """
    
    if b == 0: return None
    
    A = a
    res = 0
    pot = 0
    while abs(res*b-A)>e:
        if b>a:
            pot += 1
            a*=10
        else:
            temp = 0
            while a>b:
                a-=b
                temp += 1
            res += temp*10**-pot
    return res    
    