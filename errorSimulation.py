import numpy as np

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