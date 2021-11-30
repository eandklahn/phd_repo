import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def gauss(x,A,s,m):
    
    return A*np.exp(-(x-m)**2/(2*s**2))

def double_well(x,A,s,m):

    return -(gauss(x,A,s,m)+gauss(x,A,s,-m))

def double_well_translated(x,A,s,m,c):
    
    return c+double_well(x,A,s,m)

A = 100
s = 1.0
m = 2.6    

x = np.linspace(-10,10,1001)
y = double_well(x,A,s,m)

vals = [(0,m),(m,10)]
ymin = 97
ymax = 7
heights = -np.linspace(np.sqrt(ymax), np.sqrt(ymin), 10)**2
points = []

f, ax = plt.subplots()
for h in heights:
    points = []
    for i,bracket in enumerate(vals):
        res = root_scalar(double_well_translated,
                        args=(A,s,m,-h),
                        bracket=bracket)
        points.append(res.root)
    ax.plot((points[0], points[1]), (h,h))
    ax.plot((-points[1], -points[0]), (h,h))

ax.plot(x,y,
        c='k',
        linewidth=1)

plt.show()