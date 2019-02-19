import numpy as np

def _to_float(val):
    
    val = [float(x) for x in val.split('/')]
    
    return val[0]/val[1]

def _construct_r(r):
    
    p = {'x': np.array([1,0,0]),
         'y': np.array([0,1,0]),
         'z': np.array([0,0,1])}
    
    if '-' in r:
        return -1*p[r[-1]]
    else:
        return p[r[-1]]
        
def _construct_t(t):
    
    if '+' in t:
        return _to_float(t[1:])
    elif '-' in t:
        return -1*_to_float(t[1:])
    else:
        return 0

def _read_coordinate(s):
    
    letters = ['x', 'y', 'z']
    
    r = ''
    t = ''
    
    N, n = len(s), 0
    while n<N:
        if s[n] in letters:
            r += s[:n+1]
            t += s[n+1:]
        n += 1
    
    return [_construct_r(r), _construct_t(t)]

def _read_symmop(symm):
    
    symm = [x.strip() for x in symm.split(',')]

    R, t = np.zeros((3,3)), np.zeros(3)
    
    for n in range(3):
        L = _read_coordinate(symm[n])
        R[n,:] = L[0][:]
        t[n] = L[1]
    
    return {'R': R, 't': t}
    

if __name__ == '__main__':
    
    tests = [ 'x, y, z',
 'x, -y, z+1/2',
 'x+1/2, y+1/2, z',
 'x+1/2, -y+1/2, z+1/2']
    
    for t in tests:
    
        S = _read_symmop(t)
        
        print(S[0])
        print(S[1])