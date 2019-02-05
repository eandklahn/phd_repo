import CifFile
import numpy as np
from calculateCHILSQRs import readCHILSQ_FRs
import crystallography
import crystallographyDatabase as CD


symmOps = {'P 1 21/c 1': [
                          {'R': np.array([[1,0,0],[0,1,0],[0,0,1]]), 't': np.array([0,0,0]), 'name': 'identity'},
                          {'R': np.array([[-1,0,0],[0,-1,0],[0,0,-1]]), 't': np.array([0,0,0]), 'name': 'inversion'},
                          {'R': np.array([[-1,0,0],[0,1,0],[0,0,-1]]), 't': np.array([0,0.5,0.5]), 'name': '21 along b'},
                          {'R': np.array([[1,0,0],[0,-1,0],[0,0,1]]), 't': np.array([0,0.5,0.5]), 'name': 'c perp b'}]
                          }
            
def _search_string(inp, find):
    """ Searches the string 'inp' to locate all occurences of 'find' """
    
    occurences = []
    val = 0
    while val != -1:
        val = data.find(s, val+1)
        occurences.append(val)
    
    return tuple(occurences[:-1])

def _show_context(inp, location):
    print(inp[location-30:location+30])

def _find_loops(data):

    n = 0
    occurences = []
    while n < len(data):
        if 'data_' in data[n]:
            occurences.append(n)
        n += 1
    
    return tuple(occurences)

def _split_val_sigma(val):
    """Takes a string of the format '#1(#2)' and returns a tuple of two floats, where
    the first is #1 and the second is #2"""
    try:
        val, sig = val.split('(')[0], val.split('(')[1].split(')')[0]
        decimal_count = len(val.split('.')[1])
        sig = '0.'+'0'*(decimal_count-len(sig))+sig
    except IndexError:
        """Because sometimes there is no error on the value,
        so asking for element 1 after splitting gives an IndexError"""
        sig = 0

    return (float(val), float(sig))

def _get_type(lbl):

    type = ''
    for c in lbl:
        try:
            int(c)
        except ValueError:
            type += c

    return type

def _calculate_Fnucl(h,k,l,structure, SG, crysStruct):

    F_hkl = 0
    H = np.array(np.transpose([h,k,l]))
    
    a_ = crysStruct.aStar
    b_ = crysStruct.bStar
    c_ = crysStruct.cStar
    alpha_ = crysStruct.alphaStar
    beta_ = crysStruct.betaStar
    gamma_ = crysStruct.gammaStar
    
    for A in structure.keys():
        atom = structure[A]
        b = CD.atomInfo[atom['element']]['b_neut']
        X = atom['X']
        U = atom['U']
        for op in symmOps[SG]:
            X_current = np.matmul(op['R'],X) + op['t']
            for i in range(len(X_current)):
                while X_current[i] < 0:
                    X_current[i] += 1
                while X_current[i] > 1:
                    X_current[i] -= 1
            
            f = b*np.exp(-2j*np.pi*np.dot(H,X_current))
            
            # Debye-Waller factor
            
            Wj = 2*np.pi**2*(U[0]*h**2*a_**2 + U[1]*k**2*b_**2 + U[2]*l**2*c_**2
                             + 2*(U[3]*k*l*b_*c_*np.cos(alpha_)
                                + U[4]*h*l*a_*c_*np.cos(beta_)
                                + U[5]*h*k*a_*b_*np.cos(gamma_)))
            
            F_hkl += f*np.exp(-Wj)
            
    return F_hkl

def readCif2Structure(filename, blockname):
    
    f = open(filename)
    data = f.readlines()
    f.close()
    
    cf = CifFile.ReadCif(filename)
    block = cf[blockname]
    structure = {}
    
    #Crystallographic information
    SG = block['_space_group_name_H-M_alt']
    CS = block['_space_group_crystal_system']
    a = float(_split_val_sigma(block['_cell_length_a'])[0])
    b = float(_split_val_sigma(block['_cell_length_b'])[0])
    c = float(_split_val_sigma(block['_cell_length_c'])[0])
    alpha = float(_split_val_sigma(block['_cell_angle_alpha'])[0])
    beta = float(_split_val_sigma(block['_cell_angle_beta'])[0])
    gamma = float(_split_val_sigma(block['_cell_angle_gamma'])[0])
    P = [a, b, c, alpha, beta, gamma]
    
    #Structural information
    lbl    = block['_atom_site_label']
    x      = block['_atom_site_fract_x']
    y      = block['_atom_site_fract_y']
    z      = block['_atom_site_fract_z']
    U_type = block['_atom_site_adp_type']
    U_iso  = block['_atom_site_U_iso_or_equiv']
    
    #Thermal vibrations
    Uatomlbl = block['_atom_site_aniso_label']
    U11      = block['_atom_site_aniso_U_11']
    U22      = block['_atom_site_aniso_U_22']
    U33      = block['_atom_site_aniso_U_33']
    U23      = block['_atom_site_aniso_U_23']
    U13      = block['_atom_site_aniso_U_13']
    U12      = block['_atom_site_aniso_U_12']
    
    for i in range(len(lbl)):
        xi, dxi = _split_val_sigma(x[i])
        yi, dyi = _split_val_sigma(y[i])
        zi, dzi = _split_val_sigma(z[i])
        structure[lbl[i]] = {'X': np.transpose(np.array([xi, yi, zi])),
                            'dX': np.transpose(np.array([dxi, dyi, dzi])),
                            'Utype': U_type[i],
                            'U': [],
                            'dU': [],
                            'N': i+1,
                            'label': lbl[i],
                            'element': _get_type(lbl[i])}
        
        if U_type[i] == 'Uiso':
            structure[lbl[i]]['U'], structure[lbl[i]]['dU'] = _split_val_sigma(U_iso[i])
        elif U_type[i] == 'Uani':
            Uani_index = Uatomlbl.index(lbl[i])
            U = [0]*6
            dU = [0]*6
            U[0], dU[0] = _split_val_sigma(U11[Uani_index])
            U[1], dU[1] = _split_val_sigma(U22[Uani_index])
            U[2], dU[2] = _split_val_sigma(U33[Uani_index])
            U[3], dU[3] = _split_val_sigma(U23[Uani_index])
            U[4], dU[4] = _split_val_sigma(U13[Uani_index])
            U[5], dU[5] = _split_val_sigma(U12[Uani_index])
            structure[lbl[i]]['U'] = U
            structure[lbl[i]]['dU'] = dU                     
            
    return structure, SG, P