import numpy as np
import pandas as pd

#def lennard_compare(d1,d2):

    

def dict_compare(d1, d2):
    
    list_of_values = []
    data1 = {(h,k,l):[Fo, Fs] for h,k,l,Fo,Fs in d1}
    data2 = {(h,k,l):[Fo,Fs]+data1[(h,k,l)] for h,k,l,Fo,Fs in d2 if (h,k,l) in data1.keys()}
    list_of_values = [data2[key] for key in data2.keys()]
    
    return np.array(list_of_values)

def calc_R1(data):
    
    R1 = np.sum(np.abs(data[0]-data[2])) / np.sum(data[1])
    
    return R1
    
def compare_hkl(f1, f2):

    d1 = np.genfromtxt(fname=f1, usecols=(0,1,2,3,4))
    d2 = np.genfromtxt(fname=f2, usecols=(0,1,2,3,4))
    
    identical = []
    vals_1 = []
    vals_2 = []
    for i in d1:
        for j in d2:
            if np.array_equal(i[:3], j[:3]):
                identical.append(i[:3])
                vals_1.append(i[3:])
                vals_2.append(j[3:])
                continue
                print('Found one!')
    return identical, vals_1, vals_2

def pandas_compare(d1, d2):

    d1_df = pd.DataFrame(data=d1[:,3:], index=tuple(d1[:,:3]))
    d2_df = pd.DataFrame(data=d2[:,3:], index=tuple(d2[:,:3]))    
    
if __name__ == '__main__':
    
    f1 = '01.raw'
    f2 = '02.raw'
    
    names = ['h', 'k', 'l', 'Fo', 'Fc']
    
    d1 = np.genfromtxt(fname=f1, usecols=(0,1,2,3,4))
    d2 = np.genfromtxt(fname=f2, usecols=(0,1,2,3,4))
    
    print('Done reading')
    
    data = dict_compare(d1,d2)
    
    print(calc_R1(data))