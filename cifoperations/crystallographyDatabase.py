import numpy as np

"""
HERE WILL FOLLOW AN EXPLANATION OF WHAT IS IN THE MODULE
"""

spaceGroups = {

               'P -1'      : [
                              {'R': np.array([[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]), 't': np.array([0  ,0  , 0  ])},
                              {'R': np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]]), 't': np.array([0  ,0  , 0  ])},
                             ]                                          
                             ,
               'C 2/c'     : [
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
               'P 42/n a'  : [                                          
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

"""
b_neut: neutron scattering length for the given nucleus. Unit: 10**-14 m
"""
atomInfo = {'H' : {'atom_no': 1,  'b_neut': -0.374},
            'C' : {'atom_no': 6,  'b_neut': 0.665},
            'N' : {'atom_no': 7,  'b_neut': 0.936},
            'O' : {'atom_no': 8,  'b_neut': 0.5805},
            'S' : {'atom_no': 16, 'b_neut': 0.2847},
            'Cl': {'atom_no': 17, 'b_neut': 0.95792},
            'Co': {'atom_no': 27, 'b_neut': 0.249},
            'Dy': {'atom_no': 66, 'b_neut': 1.69}
            }