from cifoperations import crystallographyClasses as cc
import numpy as np

def UB_from_ubfrom(string):

    UB = np.array([float(s) for s in string.split()]).reshape((3,3))
    
    return UB

def mag_field_dir(cs, UB):
    
    UB_inv = np.linalg.inv(UB)
    H = np.matmul(cs.G_,
                  np.matmul(UB_inv,
                            np.array([[0],[0],[1]])
                            )
                  )
                  
    return H

dyop = cc.crystalStructure('mo_AT_467_slow_100K_0m.cif',
                           blockname='mo_at_467_slow_100k_0m')
                           
H_phi = np.array([[0],[0],[1]])

UB_1 = np.array([[-0.0143, -0.0103,  0.0450],
                 [ 0.0230, -0.0491, -0.0032],
                 [ 0.0673,  0.0146,  0.0210]])
                 
UB_2 = np.array([[-0.0567,  0.0324, -0.0104],
                 [ 0.0449,  0.0410,  0.0079],
                 [ 0.0052, -0.0003, -0.0480]])
                 
UB_3 = np.array([[ 0.0706, -0.0107,  0.0140],
                 [-0.0143, -0.0511, -0.0056],
                 [ 0.0079,  0.0029, -0.0474]])
                 
UB_3_new = np.array([[-0.0244, -0.0492, -0.0026],
                     [-0.0681,  0.0175, -0.0124],
                     [ 0.0045, -0.0018, -0.0481]])

UB_glued1 = np.array([[-0.0358,-0.0437, 0.0067],
                      [ 0.0197,-0.0211,-0.0400],
                      [ 0.0599,-0.0192, 0.0288]])
                      
UB_glued2 = np.array([[-0.0521, 0.0361,-0.0111],
                      [-0.0505,-0.0373,-0.0034],
                      [-0.0011, 0.0056, 0.0484]])
                      
UB_glued3 = np.array([[0.0281,-0.0145,-0.0392],
                      [0.0668, 0.0088, 0.0260],
                      [0.0036,-0.0494, 0.0162]])

UB_O1 = UB_from_ubfrom(''' -0.0142 -0.0103  0.0450
  0.0235 -0.0490 -0.0031
  0.0671  0.0150  0.0209
''')

UB_O21 = UB_from_ubfrom('''-0.0571  0.0321 -0.0106
  0.0444  0.0412  0.0081
  0.0056 -0.0001 -0.0479
''')

UB_O3 = UB_from_ubfrom(''' -0.0282 -0.0139  0.0394
 -0.0667  0.0085 -0.0262
 -0.0035 -0.0496 -0.0155
''')

print(mag_field_dir(dyop, UB_O1))
print(mag_field_dir(dyop, UB_O21))
print(mag_field_dir(dyop, UB_O3))

#print(mag_field_dir(dyop, UB_1))
#print(mag_field_dir(dyop, UB_2))
#print(mag_field_dir(dyop, UB_3))
#print(mag_field_dir(dyop, UB_3_new))
#print(mag_field_dir(dyop, UB_glued1))
#print(mag_field_dir(dyop, UB_glued2))
#print(mag_field_dir(dyop, UB_glued3))
