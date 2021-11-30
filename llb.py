import numpy as np

def UB_from_raf(file_path):
    
    with open(file_path, 'r') as f:
       UB = f.readlines()[1:4]
    
    UB = [l.strip().split() for l in UB]
    
    UB = [[float(val) for val in l] for l in UB]
    
    return np.array(UB)
    
def mag_field_dir(cs, UB):
    
    UB_inv = np.linalg.inv(UB)
    H = np.matmul(cs.G_,
                  np.matmul(UB_inv,
                            np.array([[0],[0],[1]])
                            )
                  )
                  
    return H