import numpy as np
from slice_image import load_image_from_scan, mask_circle2, create_circle_shape, mask_circle, add_padding_w_circles
import matplotlib.pyplot as plt
from scipy.stats import moment
import os

def sb_algorithm(image, k):
    
    # Casting image into a masked array for further processing
    masked_image = np.ma.masked_array(image, mask=np.zeros(image.shape))
    
    while True:
        # Calculate number of remaining pixels in image
        N = masked_image.size-np.sum(masked_image.mask)
        
        # Calculate background statistics
        bkgd_mean = masked_image.mean()
        bkgd_std = masked_image.std()
        bkgd_var = masked_image.var()
        bkgd_mom3 = moment(masked_image.flatten(), 3)
        bkgd_mom4 = moment(masked_image.flatten(), 4)
        s_of_variance = np.sqrt(1/N*(bkgd_mom4-(N-3)/(N-1)*bkgd_mom3))
        
        # Remove maximum pixel from background
        if abs(bkgd_var-bkgd_mean)>k*s_of_variance:
            max_index = masked_image.argmax()
            max_index = np.unravel_index(max_index, masked_image.shape)
            masked_image.mask[max_index]=1
        else:
            break
    
    return masked_image


exp = 715
scan = 73
refln_number = 3
k = 1

scanfile = 'HB3A_exp{0:04d}_scan{1:04d}.dat'.format(exp, scan)
file_beginning = 'HB3A_exp{0}_scan{1:04d}_'.format(exp, scan)
file_directory = ('C:\\Users\\Emil\\Documents\\Uddannelse\\'
                  + 'PhD\\PND susceptibility\\Co_PND_HFIR\\'
                  + 'HFIR, E18\\Data download\\HB3A\\exp715\\Datafiles\\'
                  )
files = [file_directory+f for f in os.listdir(file_directory)
         if f.startswith(file_beginning)
         ]
         
files = [(files[int(2*n)], files[int(2*n+1)])
         for n in range(int(len(files)/2))
         ]

fill_shape = create_circle_shape(radius=5)
         
image_up = load_image_from_scan(files[refln_number][0])[0]
image_dw = load_image_from_scan(files[refln_number][1])[0]

masked_up = sb_algorithm(image_up, k)
masked_dw = sb_algorithm(image_dw, k)

masked_up.mask = add_padding_w_circles(masked_up.mask, fill_shape)
            



f, ax = plt.subplots(ncols=4)
ax[0].imshow(masked_up.data)
ax[1].imshow(masked_up.mask)
ax[2].imshow(masked_dw.data)
ax[3].imshow(masked_dw.mask)

plt.show()
