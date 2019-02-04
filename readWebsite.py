import urllib.request as rq
from time import localtime, strftime

NIST_website = 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt'    
constants = {'Avogadro constant': {'tag': 'Na'},
             'Boltzmann constant': {'tag': 'kB'},
             'elementary charge': {'tag': 'e_charge'},
             'electron mass': {'tag': 'm_e'},
             'Bohr radius': {'tag': 'a0'},
             'Bohr magneton': {'tag': 'muB'},
             'atomic mass constant': {'tag': 'amu'},
             'Planck constant': {'tag': 'h'},
             'Planck constant over 2 pi': {'tag': 'hBar'},
             'speed of light in vacuum': {'tag': 'c'},
             'molar gas constant': {'tag': 'R'},
             'mag. constant': {'tag': 'mu0'},
             'electric constant': {'tag': 'epsilon_0'}}

def readTextfromSite(website):
    """
    Reads the text from the given website
    
    Input
    website: the URL of the site that should be read
    
    Output
    data: A bytearray of the site content
    """
    
    print('Opening website {}'.format(website))
    siteConnection = rq.urlopen(website)
    print('Connection open')
    data = siteConnection.read()
    siteConnection.close()
    print('Connection closed')
    
    return data

def transformBytearray(data):
    """
    Takes the bytearray returned by readTextfromSite and turns it
    into a list with each element being a line from the site
    
    Input
    data: the bytearray returned by readTextfromSite
    
    Output
    data_as_list: the content of data as a list
    """
    
    data_as_string = str(data)
    data_as_list = data_as_string.split('\\n')

    return data_as_list

def readNISTquantity(line):
    """
    Takes a line in the NIST table for scientific constants
    and returns the 'Quantity'-tag as a string
    
    Input
    line: the line from the table which should be analyzed
    
    Output
    quantity: string giving the 'Quantity'-tag of that line
    """
    
    line = line[:60].split()
    quantity = ' '.join(line)
    
    return quantity

def readNISTvalue(line):
    """
    Takes a line in the NIST table for scientific constants
    and returns the 'Value'-tag as a float
    
    Input
    line: the line from the table which should be analyzed
    
    Output
    value: the 'Value'-tag of that line
    """    
    
    line = line[60:85].split()
    value = ''
    for l in line:
        if '...' in l:
            l = l[:-3]
        value += l
        
    value = float(value)
    
    return value

def readNISTuncertainty(line):
    """
    Takes a line in the NIST table for scientific constants
    and returns the 'Uncertainty'-tag as a float
    
    Input
    line: the line from the table which should be analyzed
    
    Output
    uncertainty: the 'Uncertainty'-tag of that line
    """    

    line = line[86:110].split()
    uncertainty = ''.join(line)
    try:
        uncertainty = float(uncertainty)
        return uncertainty
    except ValueError:
        return uncertainty

def readNISTunit(line):
    """
    Takes a line in the NIST table for scientific constants
    and returns the 'Unit'-tag as a float
    
    Input
    line: the line from the table which should be analyzed
    
    Output
    unit: the 'Unit'-tag of that line
    """    
    
    line = line[110:].split()
    unit = ' '.join(line)
    
    return unit
    
def updateConstants(dict, D):
    """
    Takes in the dictionary of constants and the list of
    values read from either a text file or the NIST website
    
    Input
    dict: the dictionary holding a predefined set of constants to be read
    D: a list of lines in the NIST table of scientific constants
    
    Output
    dict: the appended dictionary now including values
    """
        
    for line in D:
        quantity = readNISTquantity(line)
        for key in dict.keys():
            if quantity == key:
                dict[key]['quantity'] = quantity
                dict[key]['value'] = readNISTvalue(line)
                dict[key]['uncertainty'] = readNISTuncertainty(line)
                dict[key]['unit'] = readNISTunit(line)

def makeNewConstantsFile(dict, filename, website=NIST_website):
    """
    Takes in the dictionary created by updateConstants and
    prints the content to a file in a way that can be imported in
    Python
    
    Input
    dict: the dictionary returned by updateConstants
    filename: name to give to the new file with constants
    website: the website that is printed in the file description
    
    Output
    this function has no output
    """
    
    text = '''"""
Scientific constants read from the format of the NIST website
{}
and output such that they can be read in as constants in Python

This file was last updated on {}
"""

'''.format(NIST_website, strftime("%a, %d %b %Y %H:%M:%S", localtime()))
    
    f = open(filename, 'w')
    
    f.write(text)
    
    for key in sorted(list(dict.keys())):
        f.write('{} = {} # +/- {} ; {} [unit: {}]\n\n'.format(dict[key]['tag'],
                                                           dict[key]['value'],
                                                           dict[key]['uncertainty'],
                                                           dict[key]['quantity'],
                                                           dict[key]['unit']))
    
    f.close()

def updateScientificConstants(constants=constants, website=NIST_website):
    """
    Updates the file 'scientificConstants.py' with the newest
    values from the NIST website
    
    Input
    website: the website that is accessed in order to obtain the
             required values
    
    Output
    this function has no output
    """
    
    D = transformBytearray(readTextfromSite(website))
    updateConstants(constants, D)
    makeNewConstantsFile(constants, 'scientificConstants.py')

