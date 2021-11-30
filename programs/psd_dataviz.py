import matplotlib.pyplot as plt
import numpy as np
import sys

def get_file_part(fname, start, end):
    """
    Returns lines numbered start, start+1, ... end-1
    from file fname
    """
    
    assert end>start
    
    part = ['']*(end-start)

    f = open(fname, 'r')
    n = 0
    while n!=start:
        f.readline()
        n+=1
    n = 0
    
    for n in range(end-start):
        part[n] = f.readline()
    
    f.close()
    
    return part

def find_frame_size(fname):
    
    f = open(fname, 'r')
    
    start, end = 0, 0
    n, line = 0, 0
    while line!='':
        line = f.readline()
        if line==' #\n': start = n
        elif line=='-\n':
            end = n
            break
        n+=1
    
    f.close()
    
    return end-start-2
    
def find_frames(fname):
    
    f = open(fname, 'r')
    
    starts = []
    n, line = 0, 0
    while line!='':
        line = f.readline()
        if line == ' #\n':
            starts.append(n)
        n+=1
    
    f.close()

    return starts

def get_images(n, fname, fbs, fs):

    iup = get_file_part(fname, fbs[n]+2, fbs[n]+fs+2)
    idw = get_file_part(fname, fbs[n]+fs+3, fbs[n]+2*fs+3)
    
    for i in range(fs):
        iup[i] = [float(val) for val in iup[i].split()]
        idw[i] = [float(val) for val in idw[i].split()]
    
    return np.array(iup), np.array(idw)

def read_file_into_structure(fname):

    f = open(fname)
    d = f.readlines()
    f.close()
    
    D = []
    idx = 0
    for line in d:
        if line==" #\n":
            D.append([[],[],[]])
            idx = 0
            continue
        elif line=="-\n":
            idx = 2
            continue
        elif idx==0:
            D[-1][idx] = line.strip().split()
            idx = 1
            continue
        else:
            pass
        D[-1][idx].append([float(val) for val in line.strip().split()])        
    
    return D

fname = r'C:\Users\emilk\Desktop\hoop\1\psd.dat'

D = read_file_into_structure(fname)