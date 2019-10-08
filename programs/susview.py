import argparse
import subprocess
import numpy as np

from CifFile import ReadCif
from cifoperations import cc

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='susViewer',
                                     description='Convert chi-tensor to appropriate coordinate system for visualization in standard molecular viewers')
    parser.add_argument('cifname')
    parser.add_argument('-b', '--blockname')
    parser.add_argument('-s', '--scale')
    parser.add_argument('-r', '--repr')
    
    _cmd_line_args = parser.parse_args()
    
    _cifname = _cmd_line_args.cifname
    _blockname = _cmd_line_args.blockname
    _type = _cmd_line_args.repr
    _scale = float(_cmd_line_args.scale)
    if _scale is None:
        _scale = 1
    
    cf = ReadCif(_cifname)
    if _blockname is None:
        _blockname = cf.keys()[0]
    structure = cc.crystalStructure(_cifname, blockname=_blockname)
        
    oMa = np.array([[structure.a_,0,0],
                    [0,structure.b_,0],
                    [0,0,structure.c_]])
    
    for atom in cf[_blockname]['_atom_site_moment.label']:
        
        _magnetic_atom_X_index = cf[_blockname]['_atom_site_moment.label'].index(atom)
        
        _X_list = structure.atoms[structure.atomdict[atom]]._magX
        _X_tensor = np.array([[_X_list[0], _X_list[5], _X_list[4]],
                              [_X_list[5], _X_list[1], _X_list[3]],
                              [_X_list[4], _X_list[3], _X_list[2]]])
        
        _T, _V = np.linalg.eig(_X_tensor)
        
        _ellipsoid = np.zeros((3,3))
        if _type is not None and _type=='repr':
            np.fill_diagonal(_ellipsoid, _T)
        else:
            np.fill_diagonal(_ellipsoid, 1/_T**2)
        
        _X_ijk = np.matmul(_V, np.matmul(_ellipsoid, np.transpose(_V)))
        _X_abc = np.matmul(structure.ABC_Mbt_IJK, np.matmul(_X_ijk, np.transpose(structure.ABC_Mbt_IJK)))
        _X_o = np.linalg.inv(np.matmul(oMa, np.matmul(_X_abc, oMa)))*_scale
        
        _magnetic_atom_U_index = cf[_blockname]['_atom_site_aniso_label'].index(atom)
        cf[_blockname]['_atom_site_aniso_U_11'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[0,0])
        cf[_blockname]['_atom_site_aniso_U_22'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[1,1])
        cf[_blockname]['_atom_site_aniso_U_33'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[2,2])
        cf[_blockname]['_atom_site_aniso_U_23'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[1,2])
        cf[_blockname]['_atom_site_aniso_U_13'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[0,2])
        cf[_blockname]['_atom_site_aniso_U_12'][_magnetic_atom_U_index] = '{:<7.5f}'.format(_X_o[0,1])
        
    _outname = _cifname[:-4]+'_magn.cif'
    outfile = open(_outname, 'w')
    outfile.write(cf.WriteOut())
    outfile.close()
    
    subprocess.Popen(['mercury', '{}'.format(_outname)])