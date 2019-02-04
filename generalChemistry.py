from atomic_info import allAtoms as aA

def molecularMass(molecule):

    atoms = []
    count = []
    
    for n in range(len(molecule)):
        if molecule[n].isupper():
            if len(atoms) != len(count):
                count.append('1')
            atoms.append(molecule[n])
        elif molecule[n].islower():
            atoms[-1] += molecule[n]
        elif float(molecule[n]):
            try: 
                float(molecule[n-1])
                count[-1] += molecule[n]
            except ValueError:
                count.append(molecule[n])
    
    if len(atoms) != len(count):
        count.append('1')
    
    mass = 0
    for n in range(len(atoms)):
        mass += aA[atoms[n]]['Mass']*float(count[n])
    
    return mass
    