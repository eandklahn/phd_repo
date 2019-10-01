import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
import os
from scipy.ndimage.filters import gaussian_filter
from skimage.feature import canny
from scipy.ndimage import measurements

def read_ub_from_scan(filename):
    """
    Reads the UB-matrix that is printed in the header of an HB3A scan-file
    and returns it as a 3x3 Numpy array
    """
    
    f = open(filename, 'r')
    CL = '' # current line being read
    while not CL.startswith('# ubmatrix'):
        CL = f.readline()
    f.close()
    
    CL = CL.split()
    CL = [float(x) for x in CL[3:]]
    
    UB = np.array([[CL[0],CL[1],CL[2]],
                   [CL[3],CL[4],CL[5]],
                   [CL[6],CL[7],CL[8]]])
    
    return UB

def read_wvln_from_scan(filename):
    """
    Reads the wavelength that is printed in the header of an HB3A scan-file
    and returns it as a float"""
    
    CL = '' # current line being read
    with open(filename, 'r') as f:
        while not CL.startswith('# wavelength'):
            CL = f.readline()
    
    return float(CL.strip().split()[-1])

def nielsen(image,k):
    
    _shape_in = image.shape
    _work = np.ma.masked_array(image, mask=False)
    _res = np.zeros(_shape_in)
    while np.abs(np.var(_work)-np.mean(_work))>k*np.std(_work):
        idx_highest = np.unravel_index(np.argmax(_work),_shape_in)
        _work.mask[idx_highest] = True
        
    f, ax = plt.subplots(ncols=2)
    _res = np.where(_work.mask, _work.data, 0)
    ax[0].imshow(image)
    ax[1].imshow(_res)
    plt.show()

def get_peak_region(image, sigma=1):
    
    peak_region_b = np.where(gaussian_filter(image, sigma)>0,True,False)
    peak_image = np.where(peak_region_b, image, 0)
    
    return peak_image, peak_region_b
    
def get_peak_border(peak_region_b):

    border_boolean_horizontal = np.zeros(peak_region_b.shape,
                                         dtype=bool)
    border_boolean_vertical = np.zeros(peak_region_b.shape,
                                         dtype=bool)
    y_pixels, x_pixels = peak_region_b.shape
    for i in range(y_pixels):
        for j in range(x_pixels):
            if peak_region_b[i,j]:
                try:
                    if not peak_region_b[i,j-1]:
                        border_boolean_horizontal[i,j] = True
                    elif not peak_region_b[i,j+1]:
                        border_boolean_horizontal[i,j] = True
                    elif not peak_region_b[i-1,j]:
                        border_boolean_vertical[i,j] = True
                    elif not peak_region_b[i+1,j]:
                        border_boolean_vertical[i,j] = True
                except IndexError:
                    # Accessing these elements can give an
                    # an error if the peak is close to the edge
                    pass
    
    border_boolean = np.logical_or(border_boolean_horizontal,
                                  border_boolean_vertical)
    
    return border_boolean

def mask_circle(image, center, radius):
    
    y_pixels, x_pixels = image.shape
    for i in range(center[0]-radius, center[0]+radius):
        for j in range(center[1]-radius, center[1]+radius):
            try:
                if np.sqrt((center[0]-i)**2
                            +(center[1]-j)**2
                            )<radius:
                        image[i,j] = True
            except IndexError:
                # Accessing these elements can give an
                # an error if the peak is close to the edge
                # (in that case, DO NOT push it)
                pass
                
def mask_circle2(image, center, fill_shape):
    
    shape_radius = int((fill_shape.shape[0]-1)/2)
    original = image[center[0]-shape_radius:center[0]+shape_radius+1,
                     center[1]-shape_radius:center[1]+shape_radius+1]
    
    image[center[0]-shape_radius:center[0]+shape_radius+1,
          center[1]-shape_radius:center[1]+shape_radius+1] = np.logical_or(original, fill_shape)
        
def add_padding_w_circles(peak_border_b, fill_shape):
    
    shape_radius=int((fill_shape.shape[0]-1)/2)
    data_shape = peak_border_b.shape
    peak_padding = np.zeros((data_shape[0]+2*shape_radius,
                             data_shape[1]+2*shape_radius),
                             dtype=bool)
    for i in range(data_shape[0]):
        for j in range(data_shape[1]):
            if peak_border_b[i,j]:
                mask_circle2(peak_padding, (i+shape_radius,j+shape_radius), fill_shape)
    
    return peak_padding[shape_radius:data_shape[0]+shape_radius, shape_radius:data_shape[1]+shape_radius]

def add_padding_to_peak(peak_region_b, count):
    
    for n in range(count):
        peak_region_b = add_peak_padding2(peak_region_b)[0]
    
    return peak_region_b   
    
def add_peak_padding2(peak_region_b):

    _shape = peak_region_b.shape
    _padding = np.full(_shape, False)
    for i in range(1,_shape[0]-1):
        for j in range(1,_shape[1]-1):
            if any([peak_region_b[i-1,j],
                    peak_region_b[i+1,j],
                    peak_region_b[i,j-1],
                    peak_region_b[i,j+1]]):
                _padding[i,j] = True
    
    return _padding, _padding!=peak_region_b
    
def add_peak_padding(peak_region_b):

    _shape = peak_region_b.shape
    _padding = np.full(_shape, False)
    for i in range(1,_shape[0]-1):
        for j in range(1,_shape[1]-1):
            if np.any(peak_region_b[i-1:i+2,j-1:j+2]):
                _padding[i,j] = True
    
    return _padding, _padding!=peak_region_b

def define_background(padded_peak):

    background = np.copy(padded_peak)
    
    while np.sum(background!=padded_peak)<np.sum(padded_peak):
        background = add_peak_padding2(background)[0]
    
    return background!=padded_peak
    
def make_too_large_background(padded_peak_region):
    
    c_radius = 10
    fill_shape = create_circle_shape(c_radius)
    border = get_peak_border(padded_peak_region)
    background = np.zeros(padded_peak_region.shape,
                          dtype=bool)
    while np.sum(padded_peak_region)>np.sum(background):
    
        background = np.logical_and(
                        add_padding_w_circles(border, fill_shape),
                        np.logical_not(padded_peak_region))
        c_radius += 1
        fill_shape = create_circle_shape(c_radius)
    
    return background

def get_intensity_from_image(image, peak_region, bkgd_region):

    # Extract pixels from image array
    peak_data = image[peak_region]
    bkgd_data = image[bkgd_region][:len(peak_data)]
    
    # Calculate Poisson-error from data, placing
    # sqrt(3) instead of zeros in error array
    peak_data_err = np.sqrt(peak_data)
    np.place(peak_data_err, peak_data_err==0, np.sqrt(3))
    
    bkgd_data_err = np.sqrt(bkgd_data)
    np.place(bkgd_data_err, bkgd_data_err==0, np.sqrt(3))
    
    # Calculate intensitites in peak and background
    peak_tot = np.sum(peak_data)
    bkgd_tot = np.sum(bkgd_data)
    
    # Calculate error in peak and background
    peak_err = np.sqrt(np.sum(peak_data_err**2))
    bkgd_err = np.sqrt(np.sum(bkgd_data_err**2))
    
    i_total = peak_tot - bkgd_tot
    s_total = np.sqrt(peak_err**2+bkgd_err**2)
    
    return i_total, s_total

def calculate_fr_error(i_up, s_up, i_dw, s_dw):
    
    fr = i_up/i_dw
    err = np.sqrt((1/i_dw)**2*s_up**2+(i_up/i_dw**2)**2*s_dw**2)
    
    return fr, err

def load_image_from_scan(filepath):

    content = ET.parse(filepath)
    root = content.getroot()
    data = root.find('Data/Detector').text
    
    x_pixels = int(root.find('Header/Number_of_X_Pixels').text)
    y_pixels = int(root.find('Header/Number_of_Y_Pixels').text)
    
    image = np.fromstring(data,
                          sep='\t',
                          dtype=np.int32
                          ).reshape((x_pixels, y_pixels)).T
    
    metadata = {}
    metadata['h'] = int(round(float(root.find('Motor_Positions/_h').text)))
    metadata['k'] = int(round(float(root.find('Motor_Positions/_k').text)))
    metadata['l'] = int(round(float(root.find('Motor_Positions/_l').text)))
    
    return image, metadata

def get_common_peak(image_up, image_dw, sigma_in=1):

    peak_boolean_up = get_peak_region(image_up, sigma=sigma_in)[1]
    peak_boolean_dw = get_peak_region(image_dw, sigma=sigma_in)[1]
    
    peak_boolean = np.logical_or(peak_boolean_up, peak_boolean_dw)
    
    return peak_boolean

def remove_border_patches2(background, padded_peak_region):
    
    y_pixels, x_pixels = padded_peak_region.shape
    
    lw, num = measurements.label(background)
    slice_tuples = measurements.find_objects(lw)
    for t in slice_tuples:
        rows = (t[0].start, t[0].stop)
        cols = (t[1].start, t[1].stop)
        if rows[0]==0:
            padded_peak_region[rows[0]:rows[1], cols[0]:cols[1]] = False
            background[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif rows[1]==y_pixels:
            padded_peak_region[rows[0]:rows[1], cols[0]:cols[1]] = False
            background[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif cols[0]==0:
            padded_peak_region[rows[0]:rows[1], cols[0]:cols[1]] = False
            background[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif cols[1]==x_pixels:
            padded_peak_region[rows[0]:rows[1], cols[0]:cols[1]] = False
            background[rows[0]:rows[1], cols[0]:cols[1]] = False

def remove_border_patches(image):
    
    y_pixels, x_pixels = image.shape
    
    lw, num = measurements.label(image)
    slice_tuples = measurements.find_objects(lw)
    for t in slice_tuples:
        rows = (t[0].start, t[0].stop)
        cols = (t[1].start, t[1].stop)
        if rows[0]==0:
            image[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif rows[1]==y_pixels:
            image[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif cols[0]==0:
            image[rows[0]:rows[1], cols[0]:cols[1]] = False
        elif cols[1]==x_pixels:
            image[rows[0]:rows[1], cols[0]:cols[1]] = False

            
def remove_small_patches(image, fill_shape):
    
    size_of_fill_shape = np.sum(fill_shape)
    lw, num = measurements.label(image)
    slice_tuples = measurements.find_objects(lw)
    
    print(slice_tuples)
    
    for n in range(1,num):
        size_of_shape = len(lw[lw==n])
        if size_of_fill_shape == size_of_shape:
            t = slice_tuples[n-1]
            rows = (t[0].start, t[0].stop)
            cols = (t[1].start, t[1].stop)
            image[rows[0]:rows[1], cols[0]:cols[1]] = False
    
    plt.imshow(lw)
    plt.show()

def keep_only_middle_patch(image):
    
    middle_y, middle_x = (image.shape[0]/2, image.shape[1]/2)
    lw, num = measurements.label(image)
    slice_tuples = measurements.find_objects(lw)
    current_min_dist = image.shape[0]
    print(current_min_dist)
    current_min_tuple = 0
    for i, t in enumerate(slice_tuples):
        rows = (t[0].start, t[0].stop)
        cols = (t[1].start, t[1].stop)
        y_slice = (rows[1]+rows[0])/2
        x_slice = (cols[1]+cols[0])/2
        print(x_slice, y_slice)
        dist = np.sqrt((middle_x-x_slice)**2+(middle_y-y_slice)**2)
        if dist < current_min_dist:
            current_min_tuple = t
    
    rows = (current_min_tuple[0].start, current_min_tuple[0].stop)
    cols = (current_min_tuple[1].start, current_min_tuple[1].stop)
    
    image[:rows[0],:] = False
    image[rows[1]:,:] = False
    image[:,:cols[0]] = False
    image[cols[1]:,:] = False
        
def create_circle_shape(radius):
    
    fill_shape = np.zeros((2*radius+1,2*radius+1), dtype=bool)
    for i in range(fill_shape.shape[0]):
        for j in range(fill_shape.shape[1]):
            if np.sqrt((i-radius)**2+(j-radius)**2)<radius:
                fill_shape[i,j] = True
    return fill_shape

def get_real_peak_region(filepath_up, filepath_dw, sigma=3, shape_in=None):

    image_up, metadata = load_image_from_scan(filepath_up)
    image_dw = load_image_from_scan(filepath_dw)[0]
    
    peak_region_boolean = get_common_peak(image_up, image_dw, sigma_in=sigma)
    peak_border_boolean = get_peak_border(peak_region_boolean)
    peak_padding = add_padding_w_circles(peak_border_boolean, shape_in)
    padded_peak_region = np.logical_or(peak_padding, peak_region_boolean)
    padded_peak_region = punch_out_border_and_small_patches(padded_peak_region, shape_in)
    background = make_too_large_background(padded_peak_region)
    background, padded_peak_region = remove_background_border_patches(background, padded_peak_region, shape_in)    
    
    return padded_peak_region
    
def get_flipping_ratio(filepath_up, filepath_dw, sigma=3, shape_in=None):
    
    make_figure = False
    
    image_up, metadata = load_image_from_scan(filepath_up)
    image_dw = load_image_from_scan(filepath_dw)[0]
    
    peak_region_boolean = get_common_peak(image_up, image_dw, sigma_in=sigma)
    peak_border_boolean = get_peak_border(peak_region_boolean)
    peak_padding = add_padding_w_circles(peak_border_boolean, shape_in)
    padded_peak_region = np.logical_or(peak_padding, peak_region_boolean)
    padded_peak_region = punch_out_border_and_small_patches(padded_peak_region, shape_in)
    background = make_too_large_background(padded_peak_region)
    background, padded_peak_region = remove_background_border_patches(background, padded_peak_region, shape_in)
    if np.sum(padded_peak_region)<1:
        metadata['rejection'] = 'no peak'
        return 0, 0, metadata
    else:
        if make_figure:
            min_val = min(image_up.min(), image_dw.min())
            max_val = max(image_up.max(), image_dw.max())
            
            f, ax = plt.subplots(ncols=2, nrows=3)
            ax[0,0].imshow(image_up, vmin=min_val, vmax=max_val)
            ax[0,1].imshow(image_dw, vmin=min_val, vmax=max_val)
            ax[1,0].imshow(np.where(padded_peak_region, image_up, 0),
                        vmin=min_val,
                        vmax=max_val)
            ax[1,1].imshow(np.where(padded_peak_region, image_dw, 0),
                        vmin=min_val,
                        vmax=max_val)
            ax[2,0].imshow(padded_peak_region)
            ax[2,1].imshow(background)
            f.tight_layout()
            plt.show()
        
        i_up, s_up = get_intensity_from_image(image_up,
                                            padded_peak_region,
                                            background)
        i_dw, s_dw = get_intensity_from_image(image_dw,
                                            padded_peak_region,
                                            background)
        fr, s_fr = calculate_fr_error(i_up, s_up, i_dw, s_dw)
        
        if i_up<0 or i_dw<0:
            metadata['rejection'] = 'negative intensity'
            return 0, 0, metadata
        elif i_up<2*s_up or i_dw<2*s_dw:
            metadata['rejection'] = 'I/s'
            return 0, 0, metadata
        else:
            return fr, s_fr, metadata

def punch_out_border_and_small_patches(image, fill_shape):
    
    fill_shape_size = np.sum(fill_shape)
    lw, num = measurements.label(image)
    patches = [np.where(lw==n, True, False) for n in range(num+1)]
    for p in patches:
        is_on_border = [np.any(p[0,:]), np.any(p[:,0]), np.any(p[-1,:]), np.any(p[:,-1])]
        #print(np.sum(p), np.sum(fill_shape))
        if np.any(is_on_border):
            #print('Punched out a border')
            image = np.where(p, False, image)
        elif np.sum(p)-fill_shape_size<0.1*fill_shape_size:
            #print('Punched out a circle')
            image = np.where(p, False, image)
    
    return image
    
def remove_background_border_patches(background, peaks, shape_in):

    lw, num = measurements.label(background)
    patches = [np.where(lw==n, True, False) for n in range(num+1)]
    for p in patches:
        # top, left, bottom, right
        is_on_border = [np.any(p[0,:]), np.any(p[:,0]), np.any(p[-1,:]), np.any(p[:,-1])]
        edges_touched = np.sum(is_on_border)
        if edges_touched<3 and edges_touched>0:
            # Border patches do not touch three borders at once
            border = get_peak_border(p)
            punch_out = add_padding_w_circles(border, shape_in)
            background = np.where(punch_out, False, background)
            peaks = np.where(punch_out, False, peaks)
    
    return background, peaks

def plot_intensity_histogram(image, start=0, normed=False):

    bin_edges = [x-0.5 for x in range(start, image.max()+2)]
    data, bins = np.histogram(image, bins=bin_edges, density=normed)
    plt.bar((bins[:-1]+bins[1:])/2, data, width=0.95)
    plt.show()

def process_scan(folder, exp, scan, file_out, mag_field, pol, radius=10):
    """
    Processes a scan from HB3A to extract flipping ratios into a file
    
    folder: directory with data files from HB3A
    exp: experiment number
    scan: scan number
    file_out: path of file that should contain final flipping ratios
    radius: radius of circle used for adding a padding to the peak region
    mag_field: magnitude of magnetic field in Tesla
    pol: polarisation value
    
    No returned values
    """
    
    scanfile = 'HB3A_exp{0:04d}_scan{1:04d}.dat'.format(exp, scan)
    file_beginning = 'HB3A_exp{0}_scan{1:04d}_'.format(exp, scan)
    
    ub_matrix = read_ub_from_scan(folder+scanfile)
    wvln = read_wvln_from_scan(folder+scanfile)
    P = np.matmul(np.linalg.inv(ub_matrix), np.array([[0],[0],[1]]))
    P = P/np.linalg.norm(P)
    
    # Getting the shape used to punch out areas of the frame
    fill_shape = create_circle_shape(radius)
    
    # Setting up file list
    files = [f for f in os.listdir(file_directory)
        if f.startswith(file_beginning)
        ]
    
    files = [(file_directory+files[n], file_directory+files[n+1]) for n in range(0,len(files),2)]
    
    reflections = []
    for i, t in enumerate(files):
        print('Extracting image data from image {}/{}'.format(i+1, len(files)), end='\r')
        reflections.append(get_flipping_ratio(t[0], t[1], shape_in=fill_shape))        
    
    with open(file_out, 'w') as f:
        f.write('#Wavelength{:>8.4f}\n'.format(wvln))
        f.write('#Orientation')
        f.write((9*'{:>8.4f}').format(*ub_matrix.T.flatten(order='F')))
        f.write('\n')
        f.write('#Polarisation{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(P[0,0],
                                                                                 P[1,0],
                                                                                 P[2,0],
                                                                                 polarisation,
                                                                                 -polarisation))
        f.write('#Magnetic_field{:>8.4f}\n'.format(mag_field))
        for r in reflections:
            if r[2].get('rejection') is None and abs(r[0]-1)>2*r[1]:
                f.write('{:>5d}{:>5d}{:>5d}{:>10.4f}{:>10.4f}\n'.format(
                        r[2]['h'],
                        r[2]['k'],
                        r[2]['l'],
                        r[0],
                        r[1])
                        )    
if __name__ == '__main__':
        
    # Experiment variables
    exp = 715
    scan = 113 # or 73
    mag_field = 0.78
    polarisation = 0.94
    
    file_directory = '''C:\\Users\\Emil\\Documents\\Uddannelse\\PhD\\pnd_susceptibility\\CoCl2(tu)4\\HFIR, E18\\Data download\\HB3A\\exp715\\Datafiles\\'''
    save_to_file = '''C:\\Users\\Emil\\Desktop\\rfl{}.ext'''.format(scan)
    
    process_scan(file_directory,
                 exp,
                 scan,
                 save_to_file,
                 mag_field,
                 polarisation)
    
    # EVERYTHING BELOW THIS LINE IS FOR TRIAL AND ERROR. Be careful about changing.
    #scanfile = 'HB3A_exp{0:04d}_scan{1:04d}.dat'.format(exp, scan)
    #file_beginning = 'HB3A_exp{0}_scan{1:04d}_'.format(exp, scan)
    #
    #ub = read_ub_from_scan(file_directory+scanfile).T.flatten(order='F')
    #fill_shape = create_circle_shape(radius)
    #
    