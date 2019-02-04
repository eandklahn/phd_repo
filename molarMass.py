
atoms = {'h': 1.008,
         'he': 4.003,
         'li': 6.94,
         'be': 9.0122,
         'b': 10.81,
         'c': 12.010,
         'n': 14.007,
         'o': 15.999,
         'f': 18.998,
         'ne': 20.180,
         'na': 22.990,
         'mg': 24.305,
         'al': 26.982,
         'si': 28.085,
         'dy': 162.50,
         'ho': 164.93}

def isInt(number):
    try:
        int(number)
        return True
    except ValueError:
        return False
        
def calcMolarMass(molecule):
    """
    Rule: A molecular formula must always start with a letter
    """
    
    # Make sure that all letters are lower case (for searching in atoms-dictionary)
    molecule = molecule.lower()

    n = 0
    d = 0
    L = len(molecule)
    splittedMolecule = []
    while n<L:
        # Splitting up the molecular formula according to where there are atom names (atom names are 3 long at most)
        if molecule[n:n+3] in atoms.keys():
            d = 3
        elif molecule[n:n+2] in atoms.keys():
            d = 2
        else:
            d = 1
        candidate = molecule[n:n+d]
        if candidate.isalpha():
            # If the candidate for a new entry in splittedMolecule is a string
            None
        elif n+d>L-1:
            # Check if looking at molecule[n+d] will give an index error.
            # In that case, we are at the final character of the string, and searching for further
            # numbers will not have an effect
            candidate = int(candidate)
        else:
            number = candidate
            while isInt(molecule[n+d]):
                number += molecule[n+d]
                d += 1
            candidate = int(number)
        splittedMolecule.append(candidate)
        n += d
    
    L = len(splittedMolecule)
    n = 0
    mass = 0
    print(splittedMolecule)
    while n<L:
        if isInt(splittedMolecule[n+1]):
            mass += atoms[splittedMolecule[n]]*splittedMolecule[n+1]
            n += 2
        else:
            mass += atoms[splittedMolecule[n]]
            n += 1
            
    return mass

#molecule = 'dy2c64h70n2o14f30'
#print(calcMolarMass(molecule))

print(2*atoms['dy']+70*atoms['c']+70*atoms['h']+2*atoms['n']+14*atoms['o']+42*atoms['f'])

    
        