import numpy as np

"""
These functions are used to calculate error propagation when handling coordinate transformations with matrices
Author: Emil A. Klahn, eklahn@chem.au.dk

It would probably be good to document these functions further and to write a specification for them.

"""


def uAvProductErrorProp(u, v, S):
    """Calculates the error in the scalar value of the multiplication u.T*A*v, where S is the error in A"""
    u = np.matrix(u).reshape(1,3)
    v = np.matrix(v).reshape(1,3)
    rows = S.shape[0]
    cols = S.shape[1]
    SUM = 0
    for i in range(rows):
        for j in range(cols):
            SUM += (u[0,i]*v[0,j]*S[i,j])**2
    return np.sqrt(SUM)
        
def transformUncertainties(M,S):
    """Transforms the uncertainties S in A, corresponding to A being in transformed via the formula A' = M*A*M.T"""
    newS = np.empty(S.shape)
    rows = S.shape[0]
    cols = S.shape[1]
    for i in range(rows):
        for j in range(cols):
            newS[i,j] = uAvProductErrorProp(M[i,:], M[j,:], S)
    return np.mat(newS)
    
def matrixSumErrorProp(T):
    """Calculates the propagated error in the result of a matrix sum, where T is a tuple of the errors on the matrices that are being added"""
    newS = np.empty(T[0].shape)
    for i in range(T[0].shape[0]):
        for j in range(T[0].shape[1]):
            SUM = 0
            for S in T:
                SUM += S[i,j]**2
            newS[i,j] = np.sqrt(SUM)
    return newS