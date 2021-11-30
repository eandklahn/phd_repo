import numpy as np

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

def purge_string_stds(string_in):
    """
    Removes all standard errors from a string, for example from a CIF
    by searching for parentheses and not promoting them to the output.
    
    Example:
    '-0.01(11)  4.33(11)  8.09(19)  -0.73(12)  3.93(17)  -3.196(65)'
    becomes
    '-0.01  4.33  8.09  -0.73  3.93  -3.196'
    
    string_in : str
        The input string from which the content of parentheses should be removed
        
    Returns
    ----------
    string_out : str
        The same string as input, just with the parentheses removed
    """
    
    n = 0
    send_on = True
    string_out = ''
    while n<len(string_in):
        if send_on:
            if string_in[n] == '(':
                send_on=False
            else:
                string_out += string_in[n]
        else:
            if string_in[n] == ')':
                send_on=True
        n+=1
    
    return string_out
    
def string_from_symmetric_matrix(M, order='11 22 33 12 13 23',
                                    dec=4):        
    """
    Contructs a space-separated string from a symmetric matrix M
    
    M : numpy array or equivalent
        The matrix containing values to be formatted to a string
        Assuming 3 x 3 array
    
    order: str
        Order that the elements should be put in the string
        21 means row 2, column 1 (row-major order)
    
    Returns
    ----------
    string_out : str
        String with the elements of M in the specified order
    """
    
    elements = {
    '11': M[0,0],
    '22': M[1,1],
    '33': M[2,2],
    '12': M[0,1],
    '13': M[0,2],
    '23': M[1,2]}
    
    order = order.split()
    ordered_elements = [str(elements[key].round(dec))for key in order]
    
    string_out = ' '.join(ordered_elements)
    
    return string_out

def symmetric_matrix_from_string(string_in,
                                 order='11 22 33 23 31 12'):
    """
    Constructs a symmetric matrix M from a space-separated string
    containing six decimal values and possibly errors in parenteses
    
    string_in : str
        The string containing values to be put in the matrix
        Must be formatted as '<val1> <val2> <val3> <val4> <val5> <val6>'
    
    order : str
        Order of the elements in string_in. This will determine where
        the individual values from string_in go in the matrix.
        21 means row 2, column 1 (row-major order)
        
    dec : int
        Number of decimal places to round the elements of M to for printing
    
    Returns
    ----------
    M : array
        A 3x3 array containing the given string as a matrix
    """
    
    # Dictionary to tell where any input must go in the matrix
    order_dict = {'11': 0, '22': 1, '33': 2, '23': 3, '32': 3,
                  '31': 4, '13': 4, '12': 5, '21': 5}
    
    # Getting the requested order of elements
    order = order.split()
    sort_list = [order_dict[val] for val in order]
    
    # Getting the elements to sort and send to matrix
    string_in = purge_string_stds(string_in)
    v = [float(s) for s in string_in.split()]
    
    # Sorting elements by indices given from order_dict
    v_sorted = [y for x,y in sorted(zip(sort_list, v), key=lambda pair: pair[0])]
    
    M = np.array([[v_sorted[0], v_sorted[5], v_sorted[4]],
                  [v_sorted[5], v_sorted[1], v_sorted[3]],
                  [v_sorted[4], v_sorted[3], v_sorted[2]]])
    
    return M

if __name__ == '__main__':
    
    string_in = '1.5978 2.9496 9.1107 -1.5714 3.4087 -4.2563'
    order = '11 33 22 13 23 12'
    
    symmetric_matrix_from_string(string_in, order=order)
