import numpy as np
import matplotlib.pyplot as plt
import imageio

def calc_dtht_m(tht_h, tht_m):
    
    diff = 0
    diff = tht_h-tht_m
    if diff>np.pi:
        diff -= 2*np.pi
    elif diff<-np.pi:
        diff += 2*np.pi
    return np.sin(diff)

# Animation options
t_max = 20
fps = 24

# Field values
tau = 1
w_h = tau
print('{:.2f} secs to make a full turn'.format(2*np.pi/w_h))

writer = imageio.get_writer('rotation.mp4', mode='I', fps=fps)

fig, ax = plt.subplots()
N = t_max*fps
t = 0
dt = 1/fps
tht_m = 0
x_h, y_h, x_m, y_m = 1,0,1,0
tht = np.linspace(0,2*np.pi,500)
x, y = np.cos(tht), np.sin(tht)
for i in range(N):
    print('{}/{}'.format(i+1,N), end='\r')
    t = i/fps
    
    tht_h = w_h*t%(2*np.pi)
    x_h, y_h = np.cos(tht_h), np.sin(tht_h)
    x_m, y_m = np.cos(tht_m), np.sin(tht_m)
    
    # Update dtht_m for next round
    dtht_m = calc_dtht_m(tht_h,tht_m)/tau
    
    # Plot
    ax.set_ylim(-2,2)
    ax.set_xlim(-2,2)
    ax.grid(axis='both')
    ax.set_aspect('equal')
    
    ax.plot(x_h, y_h, 'bo', label='field')
    ax.plot(x_m, y_m, 'ro', label='magn.')
    ax.plot(x,y,'k-')
    ax.text(1,1,s='tht_h = {:.2f}'.format(tht_h))
    ax.text(1,1.2,s='tht_m = {:.2f}'.format(tht_m))
    ax.text(1,0.8,s='dtht_m = {:.2f}'.format(dtht_m))
    
    ax.legend()
    
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
    
    # Update tht_m
    tht_m = (tht_m + dtht_m*dt)%(2*np.pi)
    