import numpy as np
from process_ac import Xp_, Xpp_
import matplotlib.pyplot as plt
from lmfit import Model, Parameters, minimize, report_fit
import sys

def objective(params, v, data):
    
    Xp_data = data[0]
    Xpp_data = data[1]
    
    Xp_model = Xp_dataset(params, v)
    Xpp_model = Xpp_dataset(params, v)
    
    Xp_resid = Xp_data - Xp_model
    Xpp_resid = Xpp_data - Xpp_model
    
    return np.concatenate([Xp_resid, Xpp_resid])

def Xp_dataset(params, v):
    tau = params['tau']
    alpha = params['alpha']
    Xt = params['Xt']
    Xs = params['Xs']
    return Xp_(v, Xs, Xt, tau, alpha)

def Xpp_dataset(params, v):
    tau = params['tau']
    alpha = params['alpha']
    Xt = params['Xt']
    Xs = params['Xs']
    return Xpp_(v, Xs, Xt, tau, alpha)

D = np.loadtxt('ccin.dat', skiprows=2)
v = D[:,0]
data_idx = int(sys.argv[1])

res = open('fit_results.txt', 'w')
for data_idx in range(1,13):
    
    Xp_data = D[:,2*(data_idx-1)+1]
    Xpp_data = D[:,2*(data_idx)]
    data = [Xp_data, Xpp_data]
    
    tau_init = (v[np.argmax(Xpp_data)]*2*np.pi)**-1
    
    Xp_model = Model(Xp_)
    Xpp_model = Model(Xpp_)
    
    fit_params = Parameters()
    fit_params.add('Xs', value=Xp_data[-1], min=0, max=np.inf)
    fit_params.add('Xt', value=Xp_data[0], min=0, max=np.inf)
    fit_params.add('tau', value=tau_init, min=0, max=np.inf)
    fit_params.add('alpha', value=0.1, min=0, max=np.inf)
    
    result_Xp = Xp_model.fit(Xp_data, fit_params, v=v)
    result_Xpp = Xpp_model.fit(Xpp_data, fit_params, v=v)
    out = minimize(objective, fit_params, args=(v, data))
    #report_fit(out.params)
    
    f, ax = plt.subplots()
    ax.plot(v, Xp_data, 'ko-')
    ax.plot(v, Xpp_data, 'ko-')
    ax.plot(v, Xp_dataset(out.params, v), 'o-')
    ax.plot(v, Xpp_dataset(out.params, v), 'c-')
    
    ax.set_xscale('log')
    result_str = '{:>20.10e} {:>20.10e} {:>20.10e} {:>20.10e}'.format(out.params['Xs'].value,out.params['Xt'].value,out.params['tau'].value,out.params['alpha'].value)
    print(result_str)
    res.write(result_str+'\n')
    f.savefig('fits/fit_{}.png'.format(data_idx))

res.close()