"""crystallographyDatabase

Intended to contain importable variables for other modules

spaceGroups (dictionary): Container for symmetry information for various space groups.
                            NOTE: Space groups will be added as they are needed
    key (string) : space group name in CIF '_space_group_name_H-M_alt'-notation
    val (list) : elements are symmetry operators given as dictionaries
        'R' (3x3 array) : rotational part of the operator
        't' (3x1 array) : translational part of the operator
                          
atomInfo (dictionary): Container for various element-specific values
                            NOTE: Elements will be added as they are needed
                            NOTE: No current notation for isotopes
    key (string) : element symbol
    val (varies) : dictionary containing implemented elements
        'atom_no' (integer) : atomic number, number of protons in nucleus
        'b_neut' (float) : neutron scattering length for the given nucleus in 10^-14 m
"""

import numpy as np


# NOTE: After 19/2/2019, this dictionary is not needed, as the module
# "symmetry" has been added, which makes it possible to read symmetry
# information directly from the strings of the cif instead of having
# to match it from this external dictionary
spaceGroups = {

               'P -1'      : [
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]), 't': np.array([0  ,0  , 0  ])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]), 't': np.array([0  ,0  , 0  ])},
                             ]                                          
                             ,
               'C 2/c nf'     : [ # NOTE: This has not been finished, hence the ' nf' in the space group name
                              {'R': np.array([[0],[0],[0]]), 't': np.array([0])}
                             ]
                             ,
               'P 1 21/c 1': [                                          
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]), 't': np.array([0  ,0  , 0  ])},
                              {'R': np.array([[-1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]), 't': np.array([0  ,0.5, 0.5])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]), 't': np.array([0  ,0  , 0  ])},
                              {'R': np.array([[ 1, 0, 0], [ 0,-1, 0], [ 0, 0, 1]]), 't': np.array([0  ,0.5, 0.5])}
                             ]
                             ,
               'P 42/n a'  : [ # NOTE: Not used at the moment. Looking up 'P42/n' will give the space group below             
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]), 't': np.array([0, 0, 0])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0, 1]]), 't': np.array([0, 0, 0])},
                              {'R': np.array([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0, 1]]), 't': np.array([1/2, 1/2, 1/2])},
                              {'R': np.array([[ 0, 1, 0], [-1, 0, 0], [ 0, 0, 1]]), 't': np.array([1/2, 1/2, 1/2])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]), 't': np.array([1/2, 1/2, 1/2])},
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]), 't': np.array([1/2, 1/2, 1/2])},
                              {'R': np.array([[ 0, 1, 0], [-1, 0, 0], [ 0, 0,-1]]), 't': np.array([0, 0, 0])},
                              {'R': np.array([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0,-1]]), 't': np.array([0, 0, 0])}
                             ]
                             ,
               'P 42/n'    : [                                          
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]), 't': np.array([0, 0, 0])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0, 1]]), 't': np.array([1/2, 1/2, 0])},
                              {'R': np.array([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0, 1]]), 't': np.array([0, 1/2, 1/2])},
                              {'R': np.array([[ 0, 1, 0], [-1, 0, 0], [ 0, 0, 1]]), 't': np.array([1/2, 0, 1/2])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]), 't': np.array([0, 0, 0])},
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]), 't': np.array([1/2, 1/2, 0])},
                              {'R': np.array([[ 0, 1, 0], [-1, 0, 0], [ 0, 0,-1]]), 't': np.array([0, 1/2, 1/2])},
                              {'R': np.array([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0,-1]]), 't': np.array([1/2, 0, 1/2])}
                             ]
}

atomInfo = {'H' : {'atom_no': 1,  'b_neut': -0.374},
            'C' : {'atom_no': 6,  'b_neut': 0.665},
            'N' : {'atom_no': 7,  'b_neut': 0.936},
            'O' : {'atom_no': 8,  'b_neut': 0.5805},
            'P' : {'atom_no': 15, 'b_neut': 0.513},
            'S' : {'atom_no': 16, 'b_neut': 0.2847},
            'Cl': {'atom_no': 17, 'b_neut': 0.9579},
            'Cu': {'atom_no': 29, 'b_neut': 0.7718},
            'Co': {'atom_no': 27, 'b_neut': 0.249},
            'Br': {'atom_no': 35, 'b_neut': 0.679},
            'I':  {'atom_no': 53, 'b_neut': 0.528},
            'Dy': {'atom_no': 66, 'b_neut': 1.69}
            }