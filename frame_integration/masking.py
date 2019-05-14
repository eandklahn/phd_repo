import numpy as np
import matplotlib.pyplot as plt
import imageio

def mask_circle(image, center, radius):
    
    y_pixels, x_pixels = image.shape
    for i in range(x_pixels):
        for j in range(y_pixels):
            if np.sqrt((center[0]-i)**2
                      +(center[1]-j)**2
                      )<radius:
                image[i,j] = True

def mask_circle_3d(image, center, radius):
    
    y_pixels, x_pixels, z_pixels = image.shape
    for i in range(x_pixels):
        for j in range(y_pixels):
            if np.sqrt((center[0]-i)**2
                      +(center[1]-j)**2
                      )<radius:
                pass
            else:
                image[j,i,:] = 0

def detect_edge_closeness(array, center, shape_in):
    
    radius = int((shape_in.shape[0]-1)/2)
    array_shape = array.shape
    distance()
    
                
def insert_shape_into_array(array, center, shape_in):
    pass
    
                
image = np.zeros((256,256),dtype=bool)
radius = 10

center = (radius, radius)
image[100,100] = True
image[200,200] = True

padding = np.zeros(image.shape, dtype=bool)

fill_shape = np.zeros((2*radius+1,2*radius+1), dtype=bool)
for i in range(fill_shape.shape[0]):
    for j in range(fill_shape.shape[1]):
        if np.sqrt((i-center[0])**2+(j-center[1])**2)<radius:
            fill_shape[i,j] = True

for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        if image[i,j]:
            print(i,j)
            print(i-radius,i+radius+1,j-radius,j+radius+1)
            padding[i-radius:i+radius+1,j-radius:j+radius+1] = fill_shape
plt.imshow(padding)
plt.show()               

#x_pixels = 256
#y_pixels = 256
#image_center = (x_pixels/2, y_pixels/2)
#radius = 10
#
#image = np.zeros((y_pixels, x_pixels), dtype=bool)
#mask_circle(image, (7,30), 15)
#mask_circle(image, (157,200), 17)
#plt.imshow(image)
#plt.show()