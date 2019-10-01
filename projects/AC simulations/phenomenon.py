import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Animation options
Nframes = 120
framesprsecond = 24

# Arrow options
A1 = 5
A2 = 3
w1 = 2
w2 = 1

fig = plt.figure()
ax = fig.add_subplot(111)
writer = imageio.get_writer('rotation.gif', mode='I', fps=framesprsecond)

t = 0

for i in range(Nframes):
    t = i/framesprsecond
    
    # "Field"-arrow
    x1 = A1*np.cos(w1*t)
    y1 = A1*np.sin(w1*t)
    # "Magnetization"-arrow
    x1 = A2*np.cos(w2*t)
    y1 = A2*np.sin(w2*t)
    
    # Plotting an arrow
    ax.quiver(0,0,x1,y1, angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(-(A1+1), A1+1)
    ax.set_ylim(-(A1+1), A1+1)
    ax.set_aspect('equal')
    
    # Recasting into numpy array to add to writer
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
