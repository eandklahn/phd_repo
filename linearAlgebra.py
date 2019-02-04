import numpy as np

def angleBetween(v1, v2):
    """Returns the angle between v1 and v2 provided that they are in the same coordinate system"""
    
    v1 = np.array(v1)
    v2 = np.array(v2)
    try:
        assert v1.shape == v2.shape
    except AssertionError:
        v1 = v1.T
    
    cosAngle = (np.dot(v1.T,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))[0,0]
    #Checking for prixomity to 1 such that Numpy does not give runtime warning
    if round(cosAngle,10) == 1:
        cosAngle = 1
    angle = np.arccos(cosAngle)
    #Making it such that the returned angle is always the sharp angle between the "vectors"
    if angle > np.pi/2:
        angle = abs(angle-np.pi)
    return np.degrees(angle)
    
def polarVector(t,p):
    """Returns a Numpy matrix object (in fact, a column vector) with the direction specified by the polar angles theta (t) and phi (p)"""
    u = np.mat([[np.sin(t)*np.cos(p)],
                [np.sin(t)*np.sin(p)],
                [np.cos(t)]         ])
    return u
    
def _general_dot(x1, x2, G=None):
    """
    Returns the general dot product for a non-orthonormal coordinate system.
    
    Inputs
    x1: column-vector of coordinates for the first vector
    x2: column-vector of coordinates for the second vector
    G: the metric matrix of the coordinate system.
    
    Outputs
    D: the scalar (dot) product of the two vectors
    """
    
    x1 = np.array(x1)
    x2 = np.array(x2)
    
    if G is None:
        G = np.identity(x1.shape[0])
    else:
        G = np.array(G)
        
    print(G)

    
        
if __name__ == '__main__':
    
    _general_dot([[1],[0],[1]], [[2],[0],[3]])
    
    
    
    

    


