import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

"""
fig = plt.figure()
ax = fig.gca(projection='3d')

x = np.linspace(-3,3,28)
y = np.linspace(-3,3,28)

X, Y = np.meshgrid(x,y)
Z = np.exp(-X**2-Y**2)

X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
triangles = mtri.Triangulation(X, Y).triangles

mycolors = np.zeros(X.shape)
for i in range(len(X)):
    if X[i]>0 and Y[i]>0:
        mycolors[i] = 1
    else:
        continue

colors = np.mean(mycolors[triangles], axis=1)

triang = mtri.Triangulation(X, Y, triangles)

collec = ax.plot_trisurf(triang, Z, cmap=cm.coolwarm, shade=False, linewidth=0)
collec.set_array(colors)
collec.autoscale()

plt.show()
"""

thetapoints = 20
phipoints = 40

fig = plt.figure()
ax = fig.gca(projection='3d')

t = np.linspace(0, np.pi, thetapoints)
p = np.linspace(0, 2*np.pi*(1-1/phipoints), phipoints)

T, P = np.meshgrid(t, p)
T, P = T.ravel(), P.ravel()

X, Y, Z = np.sin(T)*np.cos(P), np.sin(T)*np.sin(P), np.cos(T)
X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()

collec = ax.scatter(X,Y,Z)

plt.show()


"""
# Define the number of points that there should be in theta and phi
both = 25
n, m = both, both

# Constructing linspaces for theta and phi
theta = np.linspace(0, 2 * np.pi, num=n, endpoint=False)
phi   = np.linspace(np.pi*(1/(m+1)-0.5), np.pi*0.5, num=250, endpoint=False)

# Constructing a meshgrid from theta and phi and flattening both using the ravel-method
theta, phi = np.meshgrid(theta, phi)
theta, phi = theta.ravel(), phi.ravel()

# Adding the "north pole" to the data. This has coordinates (0, pi/2)
theta = np.append(theta, 0)
phi   = np.append(phi, np.pi/2)

mesh_x, mesh_y = ((np.pi*0.5 - phi)*np.cos(theta), (np.pi*0.5 - phi)*np.sin(theta))
triangles = mtri.Triangulation(mesh_x, mesh_y).triangles
x, y, z = np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta), np.sin(phi)

# Defining a custom color scalar field
vals = np.sin(6*phi) * np.sin(3*theta)
colors = np.mean(vals[triangles], axis=1)

# Plotting
fig = plt.figure()
ax = fig.gca(projection='3d')
cmap = plt.get_cmap('Blues')
triang = mtri.Triangulation(x, y, triangles)
collec = ax.plot_trisurf(triang, z, cmap=cmap, shade=False, linewidth=0.)
collec.set_array(colors)
collec.autoscale()
plt.show()
"""