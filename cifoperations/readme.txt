Scripts in this folder are

formFactors.py: 		holds functions to calculate magnetic form factors and
				Landé g-factors. Also, functions to update tables for
				electronic, j0 and j2 coefficients.

elec/magformfactor{}table.py: 	Python-dictionaries written into single files
				that can be loaded within other scripts.

crystallographyDatabase.py: 	Python-dictionaries that can be loaded into other
				scripts. At the time of writing, only a few symmetry
				operations and neutron scattering lengts for a few
				atoms.

crystallographyClasses.py:	Classes to read a cif-file and make something useful
				by casting it into a class and writing methods for that
				class.

cifoperations.py:		Similar to crystallographyClasses but older. Used for
				testing functions before making them part of the
				crystalStructure class.
				Note that there are two crystalStructure-classes in my
				libraries and that at some point they will have to be merged.
	