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

def makeNormalMatrix(M,S):

    X11 = np.random.normal(M[0,0], S[0,0])
    X22 = np.random.normal(M[1,1], S[1,1])
    X33 = np.random.normal(M[2,2], S[2,2])
    X12 = np.random.normal(M[0,1], S[0,1])
    X13 = np.random.normal(M[0,2], S[0,2])
    X23 = np.random.normal(M[1,2], S[1,2])
    
    X = np.array([[X11,X12,X13],
                  [X12,X22,X23],
                  [X13,X23,X33]])
                  
    return X
    
def normal_matrix_from_lists(m, esd):
    """
    Returns a 3x3 Numpy array with values of m(esd)
    The 6 values in m and esd are assumed to be in the order
    11, 22, 33, 12, 13, 23
    """
    
    X11 = np.random.normal(m[0], esd[0])
    X22 = np.random.normal(m[1], esd[1])
    X33 = np.random.normal(m[2], esd[2])
    X12 = np.random.normal(m[3], esd[3])
    X13 = np.random.normal(m[4], esd[4])
    X23 = np.random.normal(m[5], esd[5])
    
    X = np.array([[X11, X12, X13],[X12, X22, X23],[X13, X23, X33]])
    
    return X
    
if __name__ == '__main__':

    pass
