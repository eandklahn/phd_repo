import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

t = np.linspace(0,20,1001)
w = 3

s = np.sin(w*t)
c = np.cos(w*t)

X = 3
Xpp = 2
p=0
Xp = np.sqrt(X**2-Xpp**2)

f, ax = plt.subplots()
f.suptitle(r'$\chi$"'+'={:5.2f}, '.format(Xpp)+r'$\phi$'+'=-{:5.2f}'.format(p)+r'$\pi$')

plt.subplots_adjust(bottom=0.25)
ax.plot(t, c, label='Field')
signal, = ax.plot(t, Xp*c+Xpp*s, label='Signal')

axXpp = plt.axes([0.1,0.1, 0.8,0.03], facecolor='y')
Xppslider = Slider(axXpp, '', 0, 1, valinit=p, valstep=0.01)

def update(val):
    p = Xppslider.val
    Xpp = p*X
    Xp = X*np.sqrt(1-p**2)
    signal.set_ydata(Xp*c+Xpp*s)
    f.suptitle(r'$\chi$"'+'={:5.2f}, '.format(Xpp)+r'$\phi$'+'=-{:5.2f}'.format(p)+r'$\pi$')
    f.canvas.draw_idle()
    
Xppslider.on_changed(update)

ax.set_ylabel('Intensity (a.u.)')
ax.set_xlabel('t')
ax.legend(loc='upper right')
plt.show()