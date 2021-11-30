import argparse
import subprocess
import numpy as np

from CifFile import ReadCif
from cifoperations import cc

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='susViewer',
                                     description='Convert chi-tensor to appropriate coordinate system for visualization in standard molecular viewers')
    parser.add_argument('cifname')
    parser.add_argument('-b',
                        '--blockname',
                        help='specific name of the block to use from the cif (if multiple structures are contained in the same cif)')
    parser.add_argument('-s',
                        '--scale',
                        default=1,
                        help='arbitrary scaling of the susceptibility ellipsoid to make it fit the molecular structure')
    parser.add_argument('-r',
                        '--repr',
                        default='magn',
                        help='whether to show the magnetization ellipsoid or the representation ellipsoid (see DOI: 10.1088/0953-8984/14/38/307)')
    parser.add_argument('-m',
                        '--mercury',
                        default=False,
                        help='if this keyword is set to true, the script will try to open Mercury with the file showing the susceptibility tensor')
    
    _cmd_line_args = parser.parse_args()
    
    _cifname = _cmd_line_args.cifname
    cf = ReadCif(_cifname)
    
    _blockname = _cmd_line_args.blockname
    _type = _cmd_line_args.repr
    _scale = _cmd_line_args.scale
    _mercury = _cmd_line_args.mercury
    structure = cc.crystalStructure(_cifname, blockname=_blockname)
    
    # The following is not obsolete. It is used to find the right block, as not everything can be
    # done using the structure as loaded from the CIF.
    if _blockname is None:
        _blockname = structure.block
    
    oMa = np.array([[structure.a_,0,0],
                    [0,structure.b_,0],
                    [0,0,structure.c_]])
    
    # As the crystalStructure-class does not yet provide a method for printing a cry-file
    # from a structure, I am using the PyCIFRW-representation of the CIF as a dictionary.
    for atom in cf[_blockname]['_atom_site_moment.label']:
        
        _magnetic_atom_X_index = cf[_blockname]['_atom_site_moment.label'].index(atom)
        
        _X_list = structure.atoms[structure.atomdict[atom]]._magX
        _X_tensor = np.array([[_X_list[0], _X_list[5], _X_list[4]],
                              [_X_list[5], _X_list[1], _X_list[3]],
                              [_X_list[4], _X_list[3], _X_list[2]]])
        
        _T, _V = np.linalg.eig(_X_tensor)
        
        _ellipsoid = np.zeros((3,3))
        if _type=='magn':
            np.fill_diagonal(_ellipsoid, 1/_T**2)
        elif _type=='repr':
            np.fill_diagonal(_ellipsoid, _T)
        
        # The mathematics of this transformation can be found on pp. 206 in Giacovazzo
        _X_ijk = np.matmul(_V, np.matmul(_ellipsoid, np.transpose(_V)))
        _X_abc = np.matmul(structure.ABC_Mbt_IJK, np.matmul(_X_ijk, np.transpose(structure.ABC_Mbt_IJK))) # "line 1"
        # Question: why am I taking the inverse of the matrix calculated on the next line?
        # Answer pr. 15/7-2020: I don't know yet!
        # Answer pr. 7/3-2021: I believe it is due to the probability density being a function of the inverse of U
        # So first, the matrix chi needs to be known in the abc-system (line 1): X_abc = abc_Mbt_ijk * X_ijk * abc_Mbt_ijk^T
        # Afterwards, the matrix X_abc needs to be transformed to the dimensionless basis described in Giacovazzo
        # Finally, we are always using the inverse of U to plot data (Gia. eq. 3.B.5)
        #_X_o = np.linalg.inv(np.matmul(oMa, np.matmul(_X_abc, oMa)))*_scale
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
    
    print(f'Writing new file as {_outname}')
    if _mercury:
        subprocess.Popen(['mercury', '{}'.format(_outname)])