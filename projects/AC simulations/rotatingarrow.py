import numpy as np
import matplotlib.pyplot as plt
import imageio
from PIL import Image

Nframes = 120
framesprsecond = 24
arrowlength = 5

fig = plt.figure()
ax = fig.add_subplot(111)
writer = imageio.get_writer('rotation.gif', mode='I', fps=framesprsecond)
Elevels = [x/2 for x in range(-15,16)[::2]]
gJ = 4/3
muB = 1.5*0.001

for i in range(Nframes):

    # Making x- and y-components of the arrow
    theta = i*2*np.pi/Nframes
    x = arrowlength*np.cos(theta)
    y = arrowlength*np.sin(theta)
    
    # Plotting an arrow
    ax.quiver(0,0,x,y, angles='xy', scale_units='xy', scale=1)
    ax.quiver((arrowlength+1),0,0,y, angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(-(arrowlength+1), arrowlength+6)
    ax.set_ylim(-(arrowlength+1), arrowlength+1)
    
    # Plotting "energy levels"
    
    for mJ in Elevels:
        levelEnergy = muB*gJ*mJ*y
        if mJ > 0:
            ax.plot([arrowlength+2, arrowlength+3],[mJ**2*levelEnergy, mJ**2*levelEnergy], 'b-')
        else:    
            ax.plot([arrowlength+4, arrowlength+5],[mJ**2*levelEnergy, mJ**2*levelEnergy], 'b-')
        
    # Recasting into numpy arrow to add to writer
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.fromstring(fig.canvas.tostring_argb(), dtype=np.uint8)
    buf.shape = (h, w, 4)
    buf = np.roll(buf, 3, axis=2)
    w, h, d = buf.shape
    
    # Adding constructed image to write and clearing axes for next image
    writer.append_data(buf)
    ax.clear()
    
    #This line was needed if we wanted to make a PIL Image, but we are happy with a numpy array
    #PILimage = Image.frombytes("RGBA", (w,h), buf.tostring())
    
print("The cycle time of the GIF-image will be {} seconds\nat {} fps".format(Nframes/framesprsecond, framesprsecond))
