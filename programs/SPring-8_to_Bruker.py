# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------------------
#  'THE BEER-WARE LICENSE' (Revision 43):
#  <lkrause@chem.au.dk> wrote this file. As long as you retain this notice you can
#  do whatever you want with this stuff. If we meet some day, and you think this
#  stuff is worth it, you can buy me a beer in return.
#  ---------------------------------------------------------------------------------
import re, os, sys, glob, time, collections, multiprocessing, argparse
from datetime import datetime
import numpy as np

# try to import fabio
try:
    import fabio
except ImportError:
    print('The fabio module is needed!')
    raise SystemExit

# try to import warnings
# to get rid of numpy/fabio FutureWarnings
try:
    import warnings
    warnings.filterwarnings('ignore', category=FutureWarning)
except ImportError:
    pass

#import tkinter
try:
    no_tk = False
    from tkinter import Tk, filedialog
except ImportError:
    print('tkinter module not found!\n> We can live with that.')
    no_tk = True

# Instantiate the parser
parser = argparse.ArgumentParser()
parser.add_argument('-d', help='Debug flag', required=False, default=False, action='store_true')
parser.add_argument('-s', help='Show parameter flag', required=False, default=False, action='store_true')
parser.add_argument('-o', help='Overwrite flag', required=False, default=False, action='store_true')
parser.add_argument('-tc', help='TTH angle correction', required=False, default=0.05, type=float, nargs=1, metavar='float')
args = vars(parser.parse_args())
_DEBUG = args['d']
_SHOWPAR = args['s']
_TTH_CORR = args['tc']
_OVERWRITE = args['o']

# precompile what we need to extract
# GENERAL
det_dimx_exp = re.compile('SIZE1\s*=\s*(\d+)\s*;')
det_dimy_exp = re.compile('SIZE2\s*=\s*(\d+)\s*;')
det_size_exp = re.compile('CCD_DETECTOR_SIZE\s*=\s*(\d+\.\d+)\s*(\d+\.\d+)\s*;')
det_beam_exp = re.compile('CCD_SPATIAL_BEAM_POSITION\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;')
det_maxv_exp = re.compile('SATURATED_VALUE\s*=\s*(\d+)\s*;')
source_w_exp = re.compile('SCAN_WAVELENGTH\s*=\s*(\d+\.\d+)\s*;')
source_a_exp = re.compile('SOURCE_AMPERAGE\s*=\s*(\d+\.\d+)\s*mA\s*;')
source_v_exp = re.compile('SOURCE_VOLTAGE\s*=\s*(\d+\.\d+)\s*GeV\s*;')
# PER FRAME
goni_pos_exp = re.compile('CRYSTAL_GONIO_VALUES\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;')
goni_det_exp = re.compile('SCAN_DET_RELZERO\s*=\s*-*\d+\.\d+\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*;')
scan_rax_exp = re.compile('ROTATION_AXIS_NAME\s*=\s*(\w+)\s*;')
scan_num_exp = re.compile('SCAN_SEQ_INFO\s*=\s*\d+\s*\d+\s*(\d+)\s*;')
scan_inf_exp = re.compile('SCAN_ROTATION\s*=\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*-*\d+\.\d+\s*;')

def determineOverflow(data_in, lower_limit):
    '''
    Generates the overflow and underflow tables for the Bruker file format
    Only used for the length of these tables
    '''
    data = np.copy(data_in)
    underflow = data[data <= lower_limit]
    overflow16 = data[data > 255]
    overflow32 = overflow16[overflow16 > 65535]
    overflow_values = [underflow.astype(np.int32).shape[0], overflow16.astype(np.uint16).shape[0], overflow32.astype(np.uint32).shape[0]]
    #overflow_tables = [underflow.astype(np.int32), overflow16.astype(np.uint16), overflow32.astype(np.uint32)]
    return overflow_values#, overflow_tables, data.astype(np.uint8)

def headerToLines(header_dict):
    '''
    Converts the multi-entried dict values to single entry strings for direct
    output to Bruker frame
    '''
    return_dict = collections.OrderedDict()
    for key, value in header_dict.items():
        #logging.info(key, type(value))
        if type(value) is str:
            return_dict[key] = value
        elif isinstance(value, datetime):
            return_dict[key] = value.strftime('%d/%m/%y %H:%M:%S')
        elif isinstance(value, collections.Iterable):
            position = 0
            if type(value[0]) is str:
                strchar = ''
            elif type(value[0]) is int:
                strchar = 'd'
            elif type(value[0]) is float:
                strchar = '.6f'
            elif isinstance(value[0], np.generic):
                if np.issubdtype(value[0], float):
                    strchar = '.6f'
                elif np.issubdtype(value[0], int):
                    strchar = 'd'
            header_string = ''
            format_string = '{:<12' + strchar + '}'
            while(len(value[position:]) >  5):
                add_string = ''
                for single_value in value[position:position + 5]:
                    add_string += format_string.format(single_value)
                header_string += ''.join([add_string] + [' ']*(72 - len(add_string)) + ['\r\n'])
                position += 5
            format_string = '{:<' + str(72//(len(value[position:]) + 1)) + strchar + '}'
            for single_value in value[position:]:
                header_string += format_string.format(single_value)
            return_dict[key] = header_string
        elif (type(value) is float):# or np.issubdtype(value, float):
            return_dict[key] = '{:<.6f}'.format(value)
        elif (type(value) is int): # or
            return_dict[key] = '{:<d}'.format(value)
        elif isinstance(value, np.generic):
            if np.issubdtype(value, float):
                return_dict[key] = '{:<.6f}'.format(value)
            elif np.issubdtype(value, int):
                return_dict[key] = '{:<d}'.format(value)
        else:
            logging.info(key, value, type(value))
    return return_dict

def defaultBrukerHeader():
    '''
    Generates a standart dict for bruker 100 sfrm frames. Generated beforehand
    to have the OrderedDict in the right order without mixing standart values
    with filled-in ones
    '''
    header = collections.OrderedDict()
    header['FORMAT']  = 100                                           # Frame Format -- 86=SAXI, 100=Bruker
    header['VERSION'] = 18                                            # Header version number
    header['HDRBLKS'] = 15                                            # Header size in 512-byte blocks
    header['TYPE']    = 'Scan frame'                                  # String indicating kind of data in the frame
    header['SITE']    = 'Aarhus X-treme'                              # Site name
    header['MODEL']   = 'Huber'                                       # Diffractometer model
    header['USER']    = 'USER'                                        # Username
    header['SAMPLE']  = ''                                            # Sample ID
    header['SETNAME'] = ''                                            # Basic data set name
    header['RUN']     = 1                                             # Run number within the data set
    header['SAMPNUM'] = 1                                             # Specimen number within the data set
    header['TITLE']   = ''                                            # User comments (8 lines)
    header['NCOUNTS'] = -9999                                         # Total frame counts, Reference detector counts
    header['NOVERFL'] = [1, 0, 0]                                     # SAXI Format: Number of overflows
                                                                      # Bruker Format: #Underflows; #16-bit overfl; #32-bit overfl
    header['MINIMUM'] = -9999                                         # Minimum pixel value
    header['MAXIMUM'] = -9999                                         # Maximum pixel value
    header['NONTIME'] = 0                                             # Number of on-time events
    header['NLATE']   = 0                                             # Number of late events for multiwire data
    header['FILENAM'] = 'unknown.sfrm'                                # (Original) frame filename
    header['CREATED'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S').split()  # Date and time of creation
    header['CUMULAT'] = 1.0                                           # Accumulated exposure time in real hours
    header['ELAPSDR'] = 1.0                                           # Requested time for this frame in seconds
    header['ELAPSDA'] = 1.0                                           # Actual time for this frame in seconds
    header['OSCILLA'] = 0                                             # Nonzero if acquired by oscillation
    header['NSTEPS']  = 1                                             # steps or oscillations in this frame
    header['RANGE']   = 1.0                                           # Magnitude of scan range in decimal degrees
    header['START']   = 0.0                                           # Starting scan angle value, decimal deg
    header['INCREME'] = 1.0                                           # Signed scan angle increment between frames
    header['NUMBER']  = 1                                             # Number of this frame in series (zero-based)
    header['NFRAMES'] = 360                                           # Number of frames in the series
    header['ANGLES']  = [0.0, 0.0, 0.0, 0.0]                          # Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    header['NOVER64'] = [0, 0, 0]                                     # Number of pixels >  64K
    header['NPIXELB'] = [1, 1]                                        # Number of bytes/pixel; Number of bytes per underflow entry
    header['NROWS']   = [1024, 1]                                     # Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
                                                                      # for each mosaic tile, X varying fastest
    header['NCOLS']   = [768, 1]                                      # Number of pixels per row; number of mosaic tiles in X; dZ/dX
                                                                      # value for each mosaic tile, X varying fastest
    header['WORDORD'] = 0                                             # Order of bytes in word; always zero (0=LSB first)
    header['LONGORD'] = 0                                             # Order of words in a longword; always zero (0=LSW first
    header['TARGET']  = 'Silver'                                      # X-ray target material)
    header['SOURCEK'] = 50.0                                          # X-ray source kV
    header['SOURCEM'] = 0.880                                         # Source milliamps
    header['FILTER']  = 'Mirror'                                      # Text describing filter/monochromator setting
    header['CELL']    = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]          # Cell constants, 2 lines  (A,B,C,Alpha,Beta,Gamma)
    header['MATRIX']  = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0] # Orientation matrix, 3 lines
    header['LOWTEMP'] = [1, -17300, -6000]                            # Low temp flag; experiment temperature*100; detector temp*100
    header['ZOOM']    = [0.0, 0.0, 1.0]                               # Image zoom Xc, Yc, Mag
    header['CENTER']  = [0.0, 0.0, 0.0, 0.0]                          # X, Y of direct beam at 2-theta = 0
    header['DISTANC'] = 5.0                                           # Sample-detector distance, cm
    header['TRAILER'] = 0                                             # Byte pointer to trailer info (unused; obsolete)
    header['COMPRES'] = 'NONE'                                        # Text describing compression method if any
    header['LINEAR']  = [1.0, 0.0]                                    # Linear scale, offset for pixel values
    header['PHD']     = [0.0, 0.0]                                    # Discriminator settings
    header['PREAMP']  = [1, 1]                                        # Preamp gain setting
    header['CORRECT'] = 'UNKNOWN'                                     # Flood correction filename
    header['WARPFIL'] = 'UNKNOWN'                                     # Spatial correction filename
    header['WAVELEN'] = [0.56086, 0.55942, 0.56381]                   # Wavelengths (average, a1, a2)
    header['MAXXY']   = [1, 1]                                        # X,Y pixel # of maximum counts
    header['AXIS']    = 2                                             # Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    header['ENDING']  = [0.0, 0.5, 0.0, 0.0]                          # Setting angles read at end of scan
    header['DETPAR']  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]                # Detector position corrections (Xc,Yc,Dist,Pitch,Roll,Yaw)
    header['LUT']     = 'lut'                                         # Recommended display lookup table
    header['DISPLIM'] = [0.0, 0.0]                                    # Recommended display contrast window settings
    header['PROGRAM'] = 'Frame Conversion (python/fabio 0.6.0)'       # Name and version of program writing frame
    header['ROTATE']  = 0                                             # Nonzero if acquired by rotation (GADDS)
    header['BITMASK'] = '$NULL'                                       # File name of active pixel mask (GADDS)
    header['OCTMASK'] = [0, 0, 0, 0, 0, 0, 0, 0]                      # Octagon mask parameters (GADDS)
    header['ESDCELL'] = [0.001, 0.001, 0.001, 0.02, 0.02, 0.02]       # Cell ESD's, 2 lines (A,B,C,Alpha,Beta,Gamma)
    header['DETTYPE'] = 'UNKNOWN'                                     # Detector type
    header['NEXP']    = [1, 0, 0, 0, 0]                               # Number exposures in this frame; CCD bias level*100,;
                                                                      # Baseline offset (usually 32); CCD orientation; Overscan Flag
    header['CCDPARM'] = [0.0, 1.0, 1.0, 1.0, 65535.00]                # CCD parameters for computing pixel ESDs
    header['CHEM']    = '?'                                           # Chemical formula
    header['MORPH']   = '?'                                           # CIFTAB string for crystal morphology
    header['CCOLOR']  = '?'                                           # CIFTAB string for crystal color
    header['CSIZE']   = '?'                                           # String w/ 3 CIFTAB sizes, density, temp
    header['DNSMET']  = '?'                                           # CIFTAB string for density method
    header['DARK']    = 'NONE'                                        # Dark current frame name
    header['AUTORNG'] = [0.0, 0.0, 0.0, 0.0, 0.0]                     # Autorange gain, time, scale, offset, full scale
    header['ZEROADJ'] = [0.0, 0.0, 0.0, 0.0]                          # Adjustments to goniometer angle zeros
    header['XTRANS']  = [0.0, 0.0, 0.0]                               # Crystal XYZ translations
    header['HKL&XY']  = [0.0, 0.0, 0.0, 0.0, 0.0]                     # HKL and pixel XY for reciprocal space (GADDS)
    header['AXES2']   = [0.0, 0.0, 0.0, 0.0]                          # Diffractometer setting linear axes (4 ea) (GADDS)
    header['ENDING2'] = [0.0, 0.0, 0.0, 0.0]                          # Actual goniometer axes @ end of frame (GADDS)
    header['FILTER2'] = [0.0, 0.0, 0.0, 1.0]                          # Monochromator 2-theta, roll (both deg)
    return header
    
def convertFrame(fname):
    path_to, frame_name = os.path.split(fname)
    basename, ext = os.path.splitext(frame_name)
    frame_stem, frame_run, frame_num = basename[:-5], int(basename[-5:-3]), int(basename[-3:])
    outName = os.path.join(path_to, '{}{:>02}_{:>04}.sfrm'.format(frame_stem, frame_run, frame_num))
    
    # check if frame already exists
    if os.path.isfile(outName) and not _OVERWRITE:
        return
    
    # read in the frame
    rFrame = fabio.open(fname)
    initial_shape = rFrame.data.shape
    
    # the frame has to be flipped and rotated by 90 degrees
    rFrame.data = np.rot90(np.flipud(rFrame.data), k=1, axes=(1, 0))
    
    # get the frame saint/ready -> 1024x1024 pixel
    # fill the 981 dimension with zeros
    # cut away from both sides of the 1043 dimension
    # its not much but still no permanent solution
    # cut/pad with zeros to be able to work with the file
    #                  --y-  --x-
    #                 ( 981, 1043)
    zFrame = np.zeros((1024, 1043), dtype=np.int32)
    offset_0 = zFrame.shape[0] - rFrame.data.shape[0]
    zFrame[offset_0//2:offset_0//2 + rFrame.data.shape[0], :] = rFrame.data
    offset_1 = zFrame.shape[1] - 1024
    rFrame.data = zFrame[:, offset_1//2:offset_1//2 + 1024]
    
    # get the info file
    infFile = os.path.splitext(fname)[0] + '.inf'
    
    # check if info file exists
    if not os.path.isfile(infFile):
        print('ERROR: Info file is missing for: {}'.format(basename))
        return
    
    # extract header information
    with open(infFile) as rFile:
        infoFile = rFile.read()
        det_dimx = int(re.search(det_dimx_exp, infoFile).groups()[0])
        det_dimy = int(re.search(det_dimy_exp, infoFile).groups()[0])
        det_sizx, det_sizy = [float(i) for i in re.search(det_size_exp, infoFile).groups()]
        det_beay, det_beax = [float(i) for i in re.search(det_beam_exp, infoFile).groups()]
        det_maxv = int(re.search(det_maxv_exp, infoFile).groups()[0])
        source_w = float(re.search(source_w_exp, infoFile).groups()[0])
        source_a = float(re.search(source_a_exp, infoFile).groups()[0])
        source_v = float(re.search(source_v_exp, infoFile).groups()[0])
        goni_omg, goni_chi, goni_phi = [float(i) for i in re.search(goni_pos_exp, infoFile).groups()]
        goni_tth, goni_dxt = [float(i) for i in re.search(goni_det_exp, infoFile).groups()]
        scan_rax = str(re.search(scan_rax_exp, infoFile).groups()[0])
        scan_num = float(re.search(scan_num_exp, infoFile).groups()[0])
        scan_sta, scan_end, scan_inc, scan_exp = [float(i) for i in re.search(scan_inf_exp, infoFile).groups()]
    
    # adjust the beam center to the flipping/rotation
    det_beax = initial_shape[0] - det_beax
    det_beay = initial_shape[1] - det_beay
    
    goni_tth = goni_tth + (goni_tth * _TTH_CORR)
    
    # Gonisetup at SPring-8 to bruker conversion:
    # omega is off by 180 degrees and runs in
    # the negative direction -> -scan_inc
    goni_omg = 180 - goni_omg
    scan_inc = -scan_inc
    
    # tth and chi run in the negative direction
    goni_chi = -goni_chi
    goni_tth = -goni_tth
    
    # it is assumed that omega is the scan axis!
    scan_sta = 180 - scan_sta
    scan_end = 180 - scan_end
    
    # cut off negative outliers/bad pixels
    rFrame.data[rFrame.data < -100] = -1
    
    # scale the data to avoid underflow tables
    nexp2 = -1 * rFrame.data.min() + 1
    rFrame.data += nexp2
    
    # generate overflow tables and parameter for header
    underflow, overflow16, overflow32 = determineOverflow(rFrame.data, nexp2)
    
    # calculate detector pixel per cm, it is normalised to a 512x512 detector format
    # PILATUS3-1M pixel size is 0.172 mm 
    pix_per_512 = round((10.0 / 0.172) * (512.0 / (sum(initial_shape) / 2.0)), 6)

    # convert SPring-8 to Bruker angles
    RAXIS2BRUKER = {'Omega':1, 'Chi':2, 'Phi':0}
    axis_start = [goni_tth, goni_omg, goni_phi, goni_chi]
    axis_start[RAXIS2BRUKER[scan_rax]] = scan_sta
    axis_end = [goni_tth, goni_omg, goni_phi, goni_chi]
    axis_end[RAXIS2BRUKER[scan_rax]] = scan_end
    
    # fill known header items
    defaultHeader = defaultBrukerHeader()
    defaultHeader['NCOLS']   = [rFrame.data.shape[1], 1]# Number of pixels per row; number of mosaic tiles in X; dZ/dX
    defaultHeader['NROWS']   = [rFrame.data.shape[0], 1]# Number of rows in frame; number of mosaic tiles in Y; dZ/dY value
    # adjust the beam center for the filling/cutting of the frame
    defaultHeader['CENTER']  = [det_beax - offset_1//2, det_beay + offset_0//2, det_beax - offset_1//2, det_beay + offset_0//2]
    defaultHeader['CCDPARM'] = [0.00, 1.00, 1.00, 1.00, det_maxv]
    defaultHeader['DETPAR']  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    defaultHeader['DETTYPE'] = 'PILATUS3-1M    {}    0.00    0    0.001    0.0      0'.format(pix_per_512)
    defaultHeader['SITE']    = ['SPring-8/BL02B1']# Site name
    defaultHeader['MODEL']   = ['Synchrotron']# Diffractometer model
    defaultHeader['TARGET']  = ['Bending Magnet']# X-ray target material)
    defaultHeader['USER']    = ['USER']# Username
    defaultHeader['SOURCEK'] = [source_v]# X-ray source kV
    defaultHeader['SOURCEM'] = [source_a]# Source milliamps
    defaultHeader['WAVELEN'] = [source_w, source_w, source_w]# Wavelengths (average, a1, a2)
    defaultHeader['FILENAM'] = [basename]
    defaultHeader['CUMULAT'] = [scan_exp]# Accumulated exposure time in real hours
    defaultHeader['ELAPSDR'] = [scan_exp]# Requested time for this frame in seconds
    defaultHeader['ELAPSDA'] = [scan_exp]# Actual time for this frame in seconds
    defaultHeader['START']   = [scan_sta]# Starting scan angle value, decimal deg
    defaultHeader['ANGLES']  = axis_start# Diffractometer setting angles, deg. (2Th, omg, phi, chi)
    defaultHeader['ENDING']  = axis_end# Setting angles read at end of scan
    defaultHeader['TYPE']    = ['Generic {} Scan'.format(scan_rax)]# String indicating kind of data in the frame
    defaultHeader['DISTANC'] = [float(goni_dxt) / 10.0]# Sample-detector distance, cm
    defaultHeader['RANGE']   = [abs(scan_inc)]# Magnitude of scan range in decimal degrees
    defaultHeader['INCREME'] = [scan_inc]# Signed scan angle increment between frames
    defaultHeader['NUMBER']  = [frame_num]# Number of this frame in series (zero-based)
    defaultHeader['NFRAMES'] = [scan_num]# Number of frames in the series
    defaultHeader['AXIS']    = [RAXIS2BRUKER[scan_rax]]# Scan axis (1=2-theta, 2=omega, 3=phi, 4=chi)
    defaultHeader['LOWTEMP'] = [1, int((-273.15 + 20.0) * 100.0), -6000]# Low temp flag; experiment temperature*100; detector temp*100
    defaultHeader['NEXP']    = [1, nexp2, 0, 0, 0]
    defaultHeader['MAXXY']   = np.array(np.where(rFrame.data == rFrame.data.max()), np.float)[:, 0]
    defaultHeader['MAXIMUM'] = [np.max(rFrame.data)]
    defaultHeader['MINIMUM'] = [np.min(rFrame.data)]
    defaultHeader['NCOUNTS'] = [rFrame.data.sum(), 0]
    defaultHeader['NOVERFL'] = [-1, overflow16, overflow32]
    defaultHeader['NOVER64'] = [rFrame.data[rFrame.data > 64000].shape[0], 0, 0]
    defaultHeader['NSTEPS']  = [1]# steps or oscillations in this frame
    defaultHeader['NPIXELB'] = [1, 1]#bytes/pixel in main image, bytes/pixel in underflow table
    defaultHeader['COMPRES'] = ['NONE']#compression scheme if any
    defaultHeader['TRAILER'] = [0]#byte pointer to trailer info
    defaultHeader['PREAMP']  = [0, 0]
    defaultHeader['LINEAR']  = [1.00, 0.00]     
    defaultHeader['PHD']     = [1.00, 0.00]
    defaultHeader['OCTMASK'] = [0,0,0,1023,1023,2046,1023,1023]
    defaultHeader['DISPLIM'] = [0.0, 100.0]# Recommended display contrast window settings
    defaultHeader['CREATED'] = datetime.fromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%d %H:%M:%S').split()# use creation time of raw data!
    
    # generate a bruker sfrm and fill in data and header
    wFrame = fabio.bruker100image.Bruker100Image(data = rFrame.data, header = headerToLines(defaultHeader))
    
    # write the header to the frame
    wFrame.header = headerToLines(defaultHeader)
    # of course, this entry must be an integer ... 
    wFrame.header['HDRBLKS'] = 15
    # we need those 8 TITLE entries
    wFrame.header['TITLE'] = ' ' * 72 * 7
    
    # write the frame
    wFrame.write(outName)

if __name__ == '__main__':
    
    # Startup
    print('Convert SPring-8 Pilatus-1M data to\nBruker sfrm format     -V2018-06-24')
    print('!Warning: use for Omega scans only!')
    
    # Print startup parameters
    if _SHOWPAR:
        print(args)
    
    # user input: path to data
    while True:
        try:
            if no_tk:
                _path_to_data = input('Enter path to data: ') or '.'
            #TK, option to open the file browser window
            else:
                _path_to_data = input('Open file browser [return] or enter\na direct path: ') or 'StartFileBrowser'
                if _path_to_data == 'StartFileBrowser':
                    root = Tk()
                    root.withdraw() #don't show the Tk window
                    _path_to_data = filedialog.askdirectory(initialdir=os.path.dirname(os.getcwd()))
            fileList = sorted(glob.glob(os.path.join(_path_to_data, '*_*.tif')))
            if len(fileList) > 0:
                break
            continue
        except EOFError:
            raise SystemExit
    
    # create a pool of workers
    with multiprocessing.Pool() as pool:
        # number of tasks for progress bar
        bar_tasks = len(fileList)
        # width of the progress bar
        bar_width = 50
        # refresh of the progress bar
        bar_refresh = 0.1
        ########## DEBUG ##########
        if _DEBUG:
            convertFrame(fileList[0])
            raise SystemExit         
        ########## DEBUG ##########
        # map_async() , chunksize=1 feeds in file by file -> makes progress bar 'fluent'
        result = pool.map_async(convertFrame, fileList, chunksize=1)
        print('{:>3}%[{:{w1}}] ({:>{w2}}/{:>{w2}})\r'.format(0.0, '#'*int(0.0/(100.0/bar_width)), 0, bar_tasks, w1=bar_width, w2=len(str(bar_tasks))), end='')
        while not result.ready():
            bar_done = bar_tasks - result._number_left
            progress = int(float(bar_done)/float(bar_tasks) * 100.0)
            print('{:>3}%[{:{w1}}] ({:>{w2}}/{:>{w2}})\r'.format(progress, '#'*int(progress/(100.0/bar_width)), bar_done, bar_tasks, w1=bar_width, w2=len(str(bar_tasks))), end='')
            time.sleep(bar_refresh)
        pool.close()
        pool.join()
        print('{:>3}%[{:^{w}}] ({:>{w2}}/{:>{w2}})'.format(100, 'frame conversion complete', bar_tasks, bar_tasks, w=bar_width, w2=len(str(bar_tasks))))
