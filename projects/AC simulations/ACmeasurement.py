import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

t = np.linspace(0,2*np.pi,1000)

fig = plt.figure()
ax = fig.add_subplot(111)

writer = imageio.get_writer('M.gif', mode='I', fps=5)
i_list = list(range(21)) + list(range(1,20))[::-1]

for i in i_list:
    Mp = 0.05*i
    Mpp = 1-Mp
    M = 0
    
    Mr = Mp*np.cos(t)
    Mi = Mpp*np.sin(t)
    M = Mr + Mi
    
    ax.plot(t,Mr,'g-',label='Mr')
    ax.plot(t,Mi,'r-',label='Mi')
    ax.plot(t,M,'b-',label='M')
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(0,2*np.pi)
    ax.legend(loc='upper right')
    
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.fromstring(fig.canvas.tostring_argb(), dtype=np.uint8)
    buf.shape = (h, w, 4)
    buf = np.roll(buf, 3, axis=2)
    
    writer.append_data(buf)
    ax.clear()