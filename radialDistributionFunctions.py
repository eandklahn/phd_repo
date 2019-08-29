from scientificConstants import a0

def c_coef(n,l,j):
    if j==0:
        return 1
    else:
        return 2*((j-1)+l+1-n)/(((j-1)+1)*((j-1)+2*l+2))*c_coef(n,l,j-1)

def R_(Z,n,l,r):
    a = (a0)*10**10
    v_rho = 0
    for j in range(n-l):
        v_rho += c_coef(n,l,j)*(r/(n*a))**j
    R_n_l = 1/r*(r/(a*n))**(l+1)*np.exp(-r/(n*a))*v_rho
    return R_n_l
    
    
def P_r(Z,n,l,r):
    
    P_r = r**2*R_(Z,n,l,r)**2
    P_r = P_r/np.linalg.norm(P_r)

    return P_r
    
if __name__ == '__main__':
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    r = np.linspace(1e-14,30,1001)
    Z = 66
    #plt.plot(r, P_r(66,1,0,r), label='1s')
    #plt.plot(r, P_r(66,2,0,r), label='2s')
    #plt.plot(r, P_r(66,2,1,r), label='2p')
    #plt.plot(r, P_r(66,3,0,r), label='3s')
    #plt.plot(r, P_r(66,3,1,r), label='3p')
    #plt.plot(r, P_r(66,3,2,r), label='3d')
    #plt.plot(r, P_r(66,4,0,r), label='4s')
    #plt.plot(r, P_r(66,4,2,r), label='4d')
    #plt.plot(r, P_r(66,4,1,r), label='4p')
    plt.plot(r, P_r(66,4,3,r), label='4f')
    plt.plot(r, P_r(66,5,0,r), label='5s')
    plt.plot(r, P_r(66,5,1,r), label='5p')
    
    #plt.plot(r, P_r(25,3,2,r), label='3d')
    plt.legend()
    plt.show()
