import CifFile
import numpy as np
import os
from cifoperations import formFactors as ff
from cifoperations import crystallographyDatabase as CD
from cifoperations.symmetry import _read_symmop
from cifoperations import mj0, mj2, mj4, mj4
from crystallography import theta_from_d

def _split_val_sigma(val):
    """Takes a string of the format '#1(#2)' and returns a tuple of two floats, where
    the first is #1 and the second is #2"""
    
    try:
        val, sig = val.split('(')[0], val.split('(')[1].split(')')[0]
        decimal_count = len(val.split('.')[1])
        sig = '0.'+'0'*(decimal_count-len(sig))+sig
    
    except IndexError:
        # Because sometimes there is no error on the value,
        # so asking for element 1 after splitting gives an IndexError
        sig = 0

    return (float(val), float(sig))
    
    
class crystalStructure:

    def __init__(self, input, blockname=None):
        """
        Initializes a crystal structure.
        cifFile: the name of the CIF-file that has been opened to harvest data (if any)
        P: cell parameters of the crystal structure
        SG: the space group of the crystal structure
        """
        
        # Initial information
        self.atoms = []
        self.atomdict = {}
        self.equivPositions = []
        self.symmetrystrings = []
        
        # Information to be set if the input is a CIF
        self.cifFile = None
        self.cifData = None
        self.block = blockname
        
        # Cell parameters in direct space
        self.P = 0
        self.a, self.da         = 1,0
        self.b, self.db         = 1,0
        self.c, self.dc         = 1,0
        self.alpha, self.dalpha = 90,0
        self.beta, self.dbeta   = 90,0
        self.gamma, self.dgamma = 90,0
        
        self.handle_input(input)
        self.calculate_matrices()
        
    def handle_input(self, input):
        
        if isinstance(input, str):
            self.cifFile = input
            self._read_CIF()
        elif isinstance(input, list):
            self.a, self.b, self.c = input[:3]
            self.alpha, self.beta, self.gamma = input[3:]
        else:
            print('Crystal structure initialized without input')
        
    def calculate_matrices(self):
        
        # Cell information in direct space
        alpha, beta, gamma = np.radians(self.alpha), np.radians(self.beta), np.radians(self.gamma)
        self.G = np.mat([[self.a**2, self.a*self.b*np.cos(gamma), self.a*self.c*np.cos(beta)],
                         [self.a*self.b*np.cos(gamma), self.b**2, self.b*self.c*np.cos(alpha)],
                         [self.a*self.c*np.cos(beta), self.b*self.c*np.cos(alpha), self.c**2]])
        self.V = np.sqrt(np.linalg.det(self.G))
        
        # Cell information in reciprocal space
        self.a_ =        self.b*self.c*np.sin(alpha)/self.V
        self.b_ =        self.a*self.c*np.sin(beta)/self.V
        self.c_ =        self.a*self.b*np.sin(gamma)/self.V        
        self.alpha_ =    (np.cos(beta)*np.cos(gamma) - np.cos(alpha))/(np.sin(beta)*np.sin(gamma))
        self.beta_  =    (np.cos(alpha)*np.cos(gamma) - np.cos(beta))/(np.sin(alpha)*np.sin(gamma))
        self.gamma_ =    (np.cos(alpha)*np.cos(beta) - np.cos(gamma))/(np.sin(alpha)*np.sin(beta))
        
        self.G_ = np.linalg.inv(self.G)
        self.V_ = 1/self.V
        
        # r, q and p as they are defined in the text in my bachelor project
        r = np.cos(beta)
        
        q = ((np.cos(gamma)-np.cos(beta)*np.cos(alpha))/(np.sin(alpha)))
        
        p = np.sqrt(1-(q**2+r**2))
        
        # Matrix for coordinate transformation from crystallographic basis to orthonormal XYZ basis
        self.XYZ_Mct_ABC = np.mat([[self.a, self.b*np.cos(gamma), self.c*np.cos(beta)                                                ],
                                   [0     , self.b*np.sin(gamma), self.c*((np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma))],
                                   [0     , 0                   , self.V/(self.a*self.b*np.sin(gamma))                               ]])
        
        # Matrix for basis transformation from orthonormal CCSL basis to crystallographic basis
        # This is the basis transformation matrix (called 'M' in Giacovazzo) from the basis ijk to the basis abc
        self.ABC_Mbt_IJK = np.mat([[self.a*p, self.a*q,                  self.a*r                 ],
                                   [0,        self.b*np.sin(alpha),      self.b*np.cos(alpha)     ],
                                   [0,        0,                         self.c                   ]])
        
        # Matrix for coordinate transformation from orthonormal CCSL basis to crystallographic basis
        self.ABC_Mct_IJK = np.linalg.inv(np.transpose(self.ABC_Mbt_IJK))
        
        # Matrix for basis transformation from crystallographic basis to orthonormal CCSL basis
        self.IJK_Mbt_ABC = np.linalg.inv(self.ABC_Mbt_IJK)
        
        # Matrix for coordinate transformation from crystallographic basis to orthonormal CCSL basis
        self.IJK_Mct_ABC = np.transpose(self.ABC_Mbt_IJK)
        
        
        # Matrix for basis transformation from crystallographic basis to orthonormal XYZ basis
        self.XYZ_Mbt_ABC = np.transpose(np.linalg.inv(self.XYZ_Mct_ABC))
        
        # B-matrix (matrix that calculates coordinate transformations from reciprocal space to cartesian axes (ijk-system))
        self.B = np.matmul(self.IJK_Mct_ABC, self.G_)
    
    def _read_CIF(self):
        """Reads a structure from a CIF-file.
        If self.block is not set, the first data block in the structure is used.
        """
        
        path, file_name = os.path.split(self.cifFile)
        if path == '':
            path = os.getcwd()
            self.cifFile = path+'\\'+self.cifFile
        orig_path = os.getcwd()
        os.chdir(path)
        cf = CifFile.ReadCif(file_name)
        os.chdir(orig_path)
        
        if self.block is None: self.block = cf.keys()[0]
        data = cf[self.block]
        self.cifData = data
        _look_for_magnetic_atoms = True
        
        # Read cell parameters
        self.a, self.da         = _split_val_sigma(data['_cell_length_a'])
        self.b, self.db         = _split_val_sigma(data['_cell_length_b'])
        self.c, self.dc         = _split_val_sigma(data['_cell_length_c'])
        self.alpha, self.dalpha = _split_val_sigma(data['_cell_angle_alpha'])
        self.beta, self.dbeta   = _split_val_sigma(data['_cell_angle_beta'])
        self.gamma, self.dgamma = _split_val_sigma(data['_cell_angle_gamma'])
        
        # Read symmetry information
        if '_space_group_symop_operation_xyz' in data.keys():
            for s in data['_space_group_symop_operation_xyz']:
                self.symmetrystrings.append(s)
                self.equivPositions.append(_read_symmop(s))
        else:
            print('No symmetry information read')
            
        # Read atoms in the CIF
        labels = data['_atom_site_label']
        x      = data['_atom_site_fract_x']
        y      = data['_atom_site_fract_y']
        z      = data['_atom_site_fract_z']
        U_type = data['_atom_site_adp_type']
        U_iso  = data['_atom_site_U_iso_or_equiv']
        Occ    = data['_atom_site_occupancy']
        
        Uatomlbl = data['_atom_site_aniso_label']
        U11      = data['_atom_site_aniso_U_11']
        U22      = data['_atom_site_aniso_U_22']
        U33      = data['_atom_site_aniso_U_33']
        U23      = data['_atom_site_aniso_U_23']
        U13      = data['_atom_site_aniso_U_13']
        U12      = data['_atom_site_aniso_U_12']
        
        
        try:
            _chi_atom_lbl = data['_atom_site_moment.label']
            _chi_atom_11 = data['_atom_site_moment.chi_11']
            _chi_atom_22 = data['_atom_site_moment.chi_22']
            _chi_atom_33 = data['_atom_site_moment.chi_33']
            _chi_atom_23 = data['_atom_site_moment.chi_23']
            _chi_atom_13 = data['_atom_site_moment.chi_13']
            _chi_atom_12 = data['_atom_site_moment.chi_12']
            _chi_atom_determination = data['_atom_site_moment.determination']
        except KeyError:
            _look_for_magnetic_atoms = False
        
        for i in range(len(labels)):
            xi, dxi = _split_val_sigma(x[i])
            yi, dyi = _split_val_sigma(y[i])
            zi, dzi = _split_val_sigma(z[i])
            
            lbl = labels[i]
            # Position X is saved as a row vector
            X  = np.array([xi, yi, zi])
            dX  = np.array([dxi, dyi, dzi])
            Utype = U_type[i]
            occ, docc = _split_val_sigma(Occ[i])
            Uiso, dUiso = _split_val_sigma(U_iso[i])
            
            # Reading the vibrational parameters of the i'th atom based on the Utype
            if Utype == 'Uani':
                Uani_index = Uatomlbl.index(lbl)
                U, dU = np.zeros(6), np.zeros(6)
                U[0], dU[0] = _split_val_sigma(U11[Uani_index])
                U[1], dU[1] = _split_val_sigma(U22[Uani_index])
                U[2], dU[2] = _split_val_sigma(U33[Uani_index])
                U[3], dU[3] = _split_val_sigma(U23[Uani_index])
                U[4], dU[4] = _split_val_sigma(U13[Uani_index])
                U[5], dU[5] = _split_val_sigma(U12[Uani_index])
            
            # Atom is initialized with the arguments
            # __init__(self, label, X, dX, Utype, Uiso, dUiso, occ, docc, U=None, dU=None, element=None)
            _current_atom = Atom(lbl, X, dX, Utype, Uiso, dUiso, occ, docc, U=U, dU=dU)
            setattr(self, lbl, _current_atom)
            self.atoms.append(_current_atom)
            self.atomdict[lbl] = i
        
        if _look_for_magnetic_atoms:
            for atom in self.atoms:
                if atom.lbl in _chi_atom_lbl:
                    atom._is_magnetic = True
                    index = _chi_atom_lbl.index(atom.lbl)
                    atom._magX, atom._dmagX = np.zeros(6), np.zeros(6)
                    atom._magX[0], atom._dmagX[0] = _split_val_sigma(_chi_atom_11[index])
                    atom._magX[1], atom._dmagX[1] = _split_val_sigma(_chi_atom_22[index])
                    atom._magX[2], atom._dmagX[2] = _split_val_sigma(_chi_atom_33[index])
                    atom._magX[3], atom._dmagX[3] = _split_val_sigma(_chi_atom_23[index])
                    atom._magX[4], atom._dmagX[4] = _split_val_sigma(_chi_atom_13[index])
                    atom._magX[5], atom._dmagX[5] = _split_val_sigma(_chi_atom_12[index])
                    atom._magX_determined_by = _chi_atom_determination[index]
    
    @property
    def params(self):
        return [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
    
    def distance_between(self, atom1, atom2):
        
        try:
            atom1_X = self.__dict__[atom1].X
            atom2_X = self.__dict__[atom2].X
            diff = atom1_X - atom2_X
            x = np.matmul(diff, np.matmul(self.G, np.transpose(diff.T)).T)
            return np.sqrt(float(x))
        except KeyError as ke:
            print('{} does not exist in the structure'.format(ke))
            return None
    
    def length_of_vector(self, v):
    
        x = np.matmul(v, np.matmul(self.G, v.T))
        
        return np.sqrt(float(x))
    
    def angle_between(self, atom1, atom2, atom3=None):
        """
        Computes the angle between the bonds
        B1: atom2 --> atom1
        B2: atom2 --> atom3
        
        If atom3 is either 'a', 'b' or 'c', then the angle returned
        will be the angle between B1 and the crystallographic axis.
        
        Returns
        ----------
        angle: the requested angle given in degrees
        """
        
        axis_dict = {'a': np.array([1,0,0]),
                     'b': np.array([0,1,0]),
                     'c': np.array([0,0,1])}
        
        if isinstance(atom1, str) and isinstance(atom2, str):
            try:
                atom1 = self.__dict__.get(atom1)
                atom2 = self.__dict__.get(atom2)
                B1 = atom1.X - atom2.X
                B1_norm = self.distance_between(atom1.lbl, atom2.lbl)
            except KeyError as ke:
                print('{} is not in the structure.'.format(ke))
                return None
            
            if atom3 in axis_dict:
                """Calculate angle between B1 and unit cell axis"""
                B2 = axis_dict.get(atom3)
                B2_norm = self.__dict__.get(atom3)
            elif self.__dict__.get(atom3) is not None:
                """Calculate angle between 3 atoms"""
                atom3 = self.__dict__.get(atom3)
                B2 = atom3.X - atom2.X
                B2_norm = self.distance_between(atom2.lbl, atom3.lbl)
            else:
                print('{} is not in the structure.'.format(atom3))
                return None
        
        elif isinstance(atom1, np.ndarray) and isinstance(atom2, np.ndarray):
            B1 = atom1
            B1_norm = self.length_of_vector(atom1)
            B2 = atom2.T
            B2_norm = self.length_of_vector(atom2)
            print(B1, B2)
        
        cos_angle = np.matmul(B1, np.matmul(self.G, B2.T).T)/(B1_norm*B2_norm)
        angle = np.degrees(np.arccos(cos_angle))
        
        return float(angle)
        
    def _write_to_xyz(self, filename=None):
    
        if filename is None:
            filename = os.path.splitext(self.cifFile)[0]+'.xyz'
        
        with open(filename, 'w') as f:
            for atom in self.atoms:
                XYZ = np.matmul(self.XYZ_Mct_ABC, atom.X)[0]
                f.write('{:6s}{:>12.6f}{:>12.6f}{:>12.6f}\n'.format(atom.element,
                                                              XYZ[0,0],
                                                              XYZ[0,1],
                                                              XYZ[0,2]))
    
    def _write_to_cry(self, _cryname):
        
        Ccards = ['C {:<10.5f}{:<10.5f}{:<10.5f}{:<10.5f}{:<10.5f}{:<10.5f}\n'.format(
                   self.a, self.b, self.c, self.alpha, self.beta, self.gamma)]
        Scards = []
        Icards = []
        Qcards = []
        Lcards = []
        Fcards = []
        Acards = []
        Tcards = []
        
        Icards.extend(['I DTYP 1\n'])
        
        Lcards.extend(['L SCAL 1.000 1.000 1.000 1.000\n',
                       'L FIX SCAL 1 1 1 1\n',
                       'L MODE 5 REFI 5 WGHT 2\n'])
        
        for s in self.symmetrystrings:
            Scards.append('S '+s+'\n')
        
        _elements_present = []
        _magnetic_atoms = []
        for atom in self.atoms:
            if atom.element not in _elements_present:
                _elements_present.append(atom.element)
            
            float_repr = '{:>10.6f}'
            
            # Make the A-card
            Acardstring = ('A {: <6s}'+3*float_repr).format(atom.lbl, *(atom.X%1))
            
            if atom.Utype == 'Uiso':
                # In this case, the U-value will fit on the A-card
                Acardstring += float_repr.format(atom.Uiso)
            else:
                # In this case, a separate T-card has to be made
                Tcardstring = ('T {: <6s}  3  '+6*float_repr).format(atom.lbl, *atom.U)
                Tcardstring += '\n'
                Tcards.append(Tcardstring)
            
            if atom.occ != 1:
                sfacstring = '{:>4s}'.format(atom.element)
                occstring = ' {}'.format(atom.occ)
                uisostring = ''
                if atom.Utype != 'Uiso':
                    uisostring = float_repr.format(atom.Uiso)
                Acardstring = Acardstring + uisostring + sfacstring + occstring + '\n'
            else:
                Acardstring += '\n'
                
            if atom._is_magnetic:
                
                _lbl = atom.lbl
                
                if atom._type_magnetic_form == 'j0':
                    _exp_approx = mj0.mj0[atom.ion]
                    Fcards.append(('F {}m 2'+'{:>8.4f}'*7+'\n').format(atom.element,
                                                                _exp_approx['a'][0],
                                                                _exp_approx['b'][0],
                                                                _exp_approx['a'][1],
                                                                _exp_approx['b'][1],
                                                                _exp_approx['a'][2],
                                                                _exp_approx['b'][2],
                                                                _exp_approx['c']))
                elif atom._type_magnetic_form == 'dipole':
                    """
                    Going into this section prints F-cards for all magnetic atoms.
                    Some of these atoms may be exactly the same, and must probably
                    then be removed from the cry-file afterwards.
                    """
                    atomL = atom._angular_L
                    atomS = atom._angular_S
                    atomJ = atom._angular_J
                    s = np.linspace(0,1.6,20)
                    f = ff.magF_(atom.ion, s, L=atomL, S=atomS, J=atomJ, type='dipole')
                    s, f = list(s), list(f)
                    for n in range(0,20,5):
                        sf_temp = zip(s[n:n+5],f[n:n+5])
                        sf_temp = [l for sub in sf_temp for l in sub]
                        Fcards.append(
                        ('F {}m 3'+'{:>5.2f}{:>10.6f}'*5+'\n').format(atom.element,*sf_temp))
                    
                Qcards = ['Q {}m FORM {}\n'.format(atom.element, _lbl),
                          'Q {} CHI 0.2 0.2 0.2 0.0 0.0 0.0\n'.format(_lbl)] + Qcards
                Lcards = ['L FIX {} CH11\n'.format(_lbl),
                          'L FIX {} CH22\n'.format(_lbl),
                          'L FIX {} CH33\n'.format(_lbl),
                          'L FIX {} CH23\n'.format(_lbl),
                          'L FIX {} CH31\n'.format(_lbl),
                          'L FIX {} CH12\n'.format(_lbl),
                          'L RELA 1 1 {0} CH11 1 {0} CH33\n'.format(_lbl),
                          'L RELA 1 1 {0} CH11 1 {0} CH22\n'.format(_lbl)] + Lcards
                
            Acards.append(Acardstring)
        
        for e in _elements_present:
            Fcards.append('F {:<9s}1 {}\n'.format(e, CD.atomInfo[e]['b_neut']))
        
        Qcards = ['Q STYP PARA\n'] + Qcards
        
        cards = Ccards+Scards+Icards+Qcards+Lcards+Fcards+Acards+Tcards
        f = open(_cryname, 'w')
        for card in cards:
            f.write(card)
        f.close()
    
    def _get_params(self):
    
        return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
    
    def _show_atom(self, lbl, position=True, vibration=True):
        """
        Prints the desired information for the atom
        """
        
        if lbl in self.atomdict.keys():
            index = self.atomdict[lbl]
            X = self.atoms[index].X
            Utype = self.atoms[index].Utype
            U = self.atoms[index].U
            print('Showing information on {}'.format(lbl))
            
            if position:
                print('The position is\n{:<10.6f}{:<10.6f}{:<10.6f}'.format(*X))
            if vibration:
                if Utype == 'Uiso':
                    print('The vibrational parameter is isotropic and has the value {:<10.6f}'.format(U))
                elif Utype == 'Uani':
                    print('The vibration is anisotropic and has the values {:<10.6f}{:<10.6f}{:<10.6f}{:<10.6f}{:<10.6f}{:<10.6f}'.format(*U))
                
        else:
            print('The requested atom cannot be found in the structure.')
            
    def _show_atoms(self, type):
        
        if type == 'all':
            for atom in self.atoms:
                self._show_atom(atom.lbl)
        else:
            for atom in self.atoms:
                if atom.element == type:
                    self._show_atom(atom.lbl)
    
    def _Q_magnitude(self,H):
        
        return 2*np.pi/self._calculate_d_spacing(H)
        
    def _calculate_tht(self, H, L):
        """
        L: wavelength in angström
        """
        
        return theta_from_d(self._calculate_d_spacing(H), L)
    
    def _calculate_d_spacing(self, H):
        """
        Calculates the d-spacing for a given reflection H in this structure
        """
        
        a = self.a
        b = self.b
        c = self.c
        A = np.radians(self.alpha)
        B = np.radians(self.beta)
        Y = np.radians(self.gamma)
        h, k, l = H[0], H[1], H[2]
        
        one_over_d_dhkl_squared = (1-np.cos(A)**2-np.cos(B)**2-np.cos(Y)**2 \
                                    + 2*np.cos(A)*np.cos(B)*np.cos(Y))**(-1)*((h/a*np.sin(A))**2 \
                                    + (k/b*np.sin(B))**2 + (l/c*np.sin(Y))**2 \
                                    + 2*k*l/(b*c)*(np.cos(B)*np.cos(Y) - np.cos(A)) \
                                    + 2*l*h/(c*a)*(np.cos(Y)*np.cos(A) - np.cos(B)) \
                                    + 2*h*k/(a*b)*(np.cos(A)*np.cos(B) - np.cos(Y)))
   
        d = np.sqrt(1/one_over_d_dhkl_squared)
        
        return d

    def _calculate_Wj(self, U, Utype, R, H):
        """
        Takes the U, Utype, rotation operator and reflection coordinates and calculates
        a Debye-Waller factor for the atom
        """
    
        if Utype == 'Uiso':
            
            d = self._calculate_d_spacing(H)
            Wj = 2*np.pi**2*U/d**2
        
        elif Utype == 'Uani':
            
            h, k, l = H[0], H[1], H[2]
            U = np.array([[U[0],U[5],U[4]],[U[5],U[1],U[3]],[U[4],U[3],U[2]]])
            Uc = np.matmul(np.matmul(R, U), np.linalg.inv(R))
            
            Wj =  2*np.pi**2*(Uc[0,0]*h**2*self.a_**2 + Uc[1,1]*k**2*self.b_**2 + Uc[2,2]*l**2*self.c_**2
                        + 2*(Uc[2,1]*k*l*self.b_*self.c_*np.cos(self.alpha_)
                            + Uc[2,0]*h*l*self.a_*self.c_*np.cos(self.beta_)
                            + Uc[1,0]*h*k*self.a_*self.b_*np.cos(self.gamma_)))
        
        return Wj
     
    def _calculate_structureFactor(self, h, k, l, type='n'):
        """
        Calculates either electronic or nuclear structure factors based on the keyword 'type'.
        
        Inputs
        h, k, l: Reciprocal space coordinates of the given reflection
        type: type of structure factor to calculate (eiter 'n' for nuclear or 'e' for electronic)
        
        Outputs
        F: a complex number which is the structure factor
        """
        
        formfactors = ''
        if type == 'n':
            formfactors = 'b_neut'
        elif type == 'e':
            formfactors = 'f_elec'
        
        H = np.array([h,k,l])
        F = 0
        
        for atom in self.atoms:
            
            b = CD.atomInfo[atom.element][formfactors]
            occ = atom.occ
            
            visited = []
            
            for symmOp in self.equivPositions:
                
                R, t = symmOp['R'], symmOp['t']
                Xc = np.matmul(R, atom.X) + t
                if tuple(Xc%1) in visited:
                    continue
                else:
                    visited.append(tuple(Xc%1))
                    Wj = self._calculate_Wj(atom.U, atom.Utype, symmOp['R'], H)
                    F += b*occ*np.exp(-2j*np.pi*np.dot(H,Xc))*np.exp(-Wj)
                
        return F
    
    def _calculate_magnetic_F(self, h, k, l, Hijk):
        """
        Calculates the magnetic structure factor as a vector in the ijk coordinate system.
        This assumes that the susceptibility tensors of atoms are given in that same coordinate system
        
        Inputs
        hkl: reciprocal lattice point for which to calculate the magnetic structure factor
        Hijk: magnetic field vector in the CCSL orthonormal coordinate system
        
        Outputs
        F: the magnetic structure factor as a vector in the ijk coordinate system with unit of 10^-14 m
        """
        
        # This is the conversion factor used to convert Bohr magnetons to a scattering length
        # Still need to figure out the correct source for this conversion factor
        _conversion_factor = 0.539/2

        F = 0
        H = np.array([h,k,l])
        
        for atom in self.atoms:
            
            if atom._is_magnetic:
                
                visited = []
                
                s = 1/(2*self._calculate_d_spacing(H))
                
                _type_formfactor = atom._type_magnetic_form
                
                f = 0
                if _type_formfactor == 'j0':
                    f = ff.magF_(atom.ion, s, L=atom._angular_L, S=atom._angular_S, J=atom._angular_J)
                elif _type_formfactor == 'dipole':
                    f = ff.magF_(atom.ion, s, L=atom._angular_L, S=atom._angular_S, J=atom._angular_J, type='dipole')
                    
                chi_ijk = np.array([[atom._magX[0], atom._magX[5], atom._magX[4]],
                                    [atom._magX[5], atom._magX[1], atom._magX[3]],
                                    [atom._magX[4], atom._magX[3], atom._magX[2]]])
                    
                for symmOp in self.equivPositions:
                    R, t = symmOp['R'], symmOp['t']
                    Xc = np.matmul(R, atom.X) + t
                    if tuple(Xc%1) in visited:
                        continue
                    else:
                        visited.append(tuple(Xc%1))
                        Rijk = np.matmul(np.linalg.inv(np.transpose(self.IJK_Mbt_ABC)),np.matmul(R,np.transpose(self.IJK_Mbt_ABC)))
                        m = np.matmul(np.matmul(np.matmul(Rijk, chi_ijk), np.linalg.inv(Rijk)), Hijk)
                        Wj = self._calculate_Wj(atom.U, atom.Utype, symmOp['R'], H)
                        F += f*m*np.exp(-2j*np.pi*np.dot(H,Xc))*np.exp(-Wj)
        
        F *= _conversion_factor
    
        return F
    
    def _calculate_Fm_perp_Q(self, h, k, l, Hijk):
        """
        Calculates the component of Fm that is perpendicular to the scattering vector
        """
        Q = np.array([[h],[k],[l]])
        
        Q_ijk = np.matmul(self.IJK_Mct_ABC,np.matmul(self.G_, Q))
        Fm_ijk = self._calculate_magnetic_F(h,k,l, Hijk)
        
        Fm_para_q_ijk_num = np.dot(np.transpose(Fm_ijk),Q_ijk)
        Fm_para_q_ijk_den = np.dot(np.transpose(Q_ijk), Q_ijk)
        
        Fm_para_q_ijk = (Fm_para_q_ijk_num/Fm_para_q_ijk_den)[0,0]*Q_ijk
        
        Fm_perp_q_ijk = Fm_ijk - Fm_para_q_ijk
        
        return Fm_perp_q_ijk
    
    def _cellParams(self):
        """
        Reports on the cell parameters of the structure
        """
        
        print('Cell parameters for this compound are\n{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'.format(
               self.a, self.b, self.c, self.alpha, self.beta, self.gamma))
               
    def _atoms(self):
        """
        Reports on all the atoms in the structure
        """
        
        for atom in self.atoms:
            print(atom.lbl)


class Atom:

    def __init__(self, label, X, dX, Utype, Uiso, dUiso, occ, docc, U=None, dU=None, element=None):
        """
        Class used to represent an atom in a crystal structure
        
        ...
        
        Attributes
        ----------
        label : string
            the label that the atom has in a structure.
        element : string
            the element of the atom
        X : array
            position of atom in unit cell in fractional coordinates 
        dX : array
            sigmas on the atom position
        Utype : string
            whether vibrations are given as Uiso or Uani
        U : array/float
            values for the vibrational parameters. Saved according to Utype
        dU : array/float
            sigmas on vibrational parameters
        
        Methods
        ----------
        _get_element
            Reads the element of the atom from the given label
        """
        
        self.lbl = label
        self.element = self._get_element()
        self.X = X
        self.dX = dX
        self.Utype = Utype
        self.Uiso = Uiso
        self.dUiso = dUiso
        self.U = U
        self.dU = dU
        self.occ = occ
        self.docc = docc
        
        # Magnetic stuff
        self._is_magnetic = False
        self.charge = 0
        self.ion = ''
        self._type_magnetic_form = 'j0'
        self._angular_L = 0
        self._angular_S = 0
        self._angular_J = 0
        
    
    def _get_element(self):
        s = self.lbl
        e = ''
        n = 0
        while True:
            if n>(len(s)-1):
                break
            elif s[n].isalpha():
                e+=s[n]
                n+=1
            else:
                break
        return e
    
    @property
    def magX_tensor(self):
        return np.array([[self._magX[0], self._magX[5], self._magX[4]],
                         [self._magX[5], self._magX[1], self._magX[3]],
                         [self._magX[4], self._magX[3], self._magX[2]]])
    
if __name__ == '__main__':    

    from calculateCHILSQRs import readCHILSQ_FRs

    filename = '20K.cif'
    block = 'phys_i13'
    crys = crystalStructure(cifFile=filename, blockname=block)
    
    data = readCHILSQ_FRs('chilsqtBu.lis')
    data = data[list(data.keys())[0]]
     
     
    for n in range(len(data['h'])):
        h            = data['h'][n]
        k            = data['k'][n]
        l            = data['l'][n]
        Fnucl_CHILSQ = data['Fmag'][n]
        Fnucl_script = crys._calculate_structureFactor(h,k,l,type='n')
        
        if Fnucl_script.imag < 1e-12:
            print('{:6d} {:6d} {:6d} {:>10.4f} {:>10.4f}'.format(int(h),int(k),int(l), Fnucl_CHILSQ, Fnucl_script.real))
        else:
            print('There is a significant imaginary contribution from reflection ({},{},{})'.format(int(h),int(k),int(l)))
    
