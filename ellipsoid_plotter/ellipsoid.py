import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_ellipsoid(X, trans=[0,0,0], axes_in=None):
    
    E, V = np.linalg.eig(X)
    
    try:
        assert (np.all(E>0))
    except AssertionError:
        print('_plot_ellipsoid closed, because not all eigenvalues are positive')
        print(E)
        return 1
    
    t = np.linspace(0,np.pi,10)
    p = np.linspace(0,2*np.pi,20)
    T, P = np.meshgrid(t, p)
    
    X = E[0]*np.sin(T)*np.cos(P)
    Y = E[1]*np.sin(T)*np.sin(P)
    Z = E[2]*np.cos(T)
    _shape = X.shape
    _maxval = np.max(E)*1.1
    
    H = np.array([X.flatten(),Y.flatten(),Z.flatten()])
    _H = np.matmul(V,H)
    _X, _Y, _Z = np.reshape(_H[0,:],_shape)+trans[0], np.reshape(_H[1,:],_shape)+trans[1], np.reshape(_H[2,:],_shape)+trans[2]
    
    if axes_in is None:
        fig = plt.figure()
        axes_in = fig.add_subplot(111, projection='3d')
        
    axes_in.plot_surface(_X,_Y,_Z,
                         facecolor=(0,1,1),
                         alpha=0.5,
                         linewidth=1
                         )
    
    # Keep these lines. They keep a reasonable aspect ratio between the axes
    axes_in.set_xlim(trans[0]-_maxval,trans[0]+_maxval)
    axes_in.set_ylim(trans[1]-_maxval,trans[1]+_maxval)
    axes_in.set_zlim(trans[2]-_maxval,trans[2]+_maxval)
    
    axes_in.set_xlabel('x')
    axes_in.set_ylabel('y')
    axes_in.set_zlabel('z')
    
    return axes_in

    
if __name__=='__main__':

    _vals = [2.64,7.07,2.44,-3.36,1.34,-3.41]
    t = [10,20,30]
    X = np.array([[_vals[0],_vals[5],_vals[4]],
                  [_vals[5],_vals[1],_vals[3]],
                  [_vals[4],_vals[3],_vals[2]]])
                  
    _plot_ellipsoid(X, t)