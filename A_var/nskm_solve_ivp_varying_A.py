#1d varying A proper
#very similar code to "odeint_nskm" file 

import numpy as np
from scipy.integrate import solve_ivp as si
import matplotlib.pyplot as plt

def A(t,Q,R,k,Am):
    return Am+ (k*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R)))

def A2(t,Q,R,k1,k2,Am):
    return Am+ (k1*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R))) + (k2*(1 + np.tanh((Q*np.sin((0.5*np.pi*t)+(np.pi))) -R)))


def ivp_nskm(t,z,Q,R,k,k2 = 0.5,Am = 0.7,B=0.45,rf = 1):
    u = z[0]
    w = z[1]
    dudt = w*(u**2) - B*u
    if rf == 1:
        dwdt = A(t,Q,R,k,Am) - w - w*(u**2)
    else:
        dwdt = A2(t,Q,R,k,k2,Am) - w - w*(u**2)
    return [dudt,dwdt]


def basic_fig_plot(z,tmax = 20,rf = 1,Q=1.3518,R=1,k=1,k2 = 0.5,Am = 0.7,B=0.45,save = False,title='plot',tmin = 0):
    '''z should be (u,w) in form tuple,list or np array''' 
    #points it will evaluate at (these are not all points it is evaluated at just the ones you want out)
    #explicitly giving it time points to get out gives a smoother graph
    tsp = np.linspace(tmin,tmax,1000)
    #solve
    sola = si(ivp_nskm,(tmin,tmax),z,t_eval = tsp,args = (Q,R,k,k2,Am,B,rf))
    #as z was a list length 2 where u = z[0] and w = z[1] can read solns like this
    u = sola.y[0]
    w = sola.y[1]
    #plot against time points
    plt.plot(tsp,u,'g')
    plt.plot(tsp,w,'b')
    plt.margins(0.01)
    plt.xlabel('time' , fontsize = 15)
    plt.ylabel('density' , fontsize = 15)
    plt.legend(['u(t)','w(t)'])
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    
    if save:
        plt.savefig('plots/' + title + '.png')
    plt.show()
    
def t_min_find(z=(0.5,0.5),tr = (0,4),tmax = 50,points=101):
    '''z is initial cond,
tr is range of values tmin can take (i.e., how far through the cycle it is
tmax is the cut off point for the solver
points is the no of evaluations'''
    tmin_vals = np.linspace(tr[0],tr[1],points)
    #
    a_val  = []
    a2_val = []
    #for each val
    for a in tmin_vals:
        #evaluate using solver
        a_sol  = si(ivp_nskm,(a,tmax+a),z,args = (1.3513,1,1,0.5,0.7,0.45,1))
        a2_sol = si(ivp_nskm,(a,tmax+a),z,args = (4,3.42,1,0.5,0.7,0.45,2))
        #if at late times is non zero  add to list of valid regions
        if a_sol.y[0][-1] > 1e-2:
            a_val.append(a)
        if a2_sol.y[0][-1] > 1e-2:
            a2_val.append(a)
    plt.plot([0,0],[0,1.5],linewidth = 0,markersize=0)#quicker than setting range
    plt.plot(a_val,np.ones(len(a_val)),linewidth = 0,marker = 's',markersize = 1)
    plt.plot(a2_val,np.ones(len(a2_val))*0.5,linewidth = 0,marker = 's',markersize = 1)
    plt.margins(0)
    plt.show()
    return a_val ,  a2_val
