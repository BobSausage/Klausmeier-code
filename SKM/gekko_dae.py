#trying to use gekko to solve the DAE

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from gekko import GEKKO
from random import randint as ri

def run_gek_oneD():
    m = GEKKO() # create model
    A = 1.1776
    B = 0.45
    u = m.Var(0.5) #create variable with initial condition
    w = m.Var(0.5) #
    m.Equation(u.dt()==(u**2)*w - (B*u))    #u equation
    m.Equation(0 == A - w - w*(u**2))       #w equation (w \dot{w} = 0)
    m.time = np.linspace(0,20,100)              #set time points

    print('initialised')
    #solve ODE
    m.options.IMODE= 4 # 1 for all derivs are 0, 4 for simultaneous solving, 7 for sequential solving
    m.solve(disp=False) # solve (do not display full output)

    print('solved')
    #plot results
    plt.plot(m.time,u,color = 'g')
    plt.plot(m.time,w,color = 'b')
    plt.show()
def run_gek_twoD_notkm():
    m = GEKKO()
    k = 1        # constant
    y = m.Array(m.Var,(4)) # create GEKKO variable
    print(y)
    for i in range(len(y)):
        y[i].value = i
    print(y)
    m.Equations([y[i].dt()==-k*y[i] for i in range(len(y))]) # create GEKKO equation
    m.time = np.linspace(0,20) # time points

    # solve ODE
    m.options.IMODE = 4
    m.solve(disp=False,debug=True)

    ##print(y)
    # plot results
    for a in range(len(y)):
        plt.plot(m.time,y[a])
    plt.xlabel('time')
    plt.ylabel('y(t)')
    plt.show()

from optimal_L import *
from L1 import *
 

###this might not be the most efficient way to do this but hopefully it works
def gek_run(n):
    #create model
    m = GEKKO(remote = False)
    #define parameters
    A = 0.9
    B = 0.45
    dis  = 0.25 #coeff of 2nd order \nabla  = 1/2h
    dis2 = 0.25  #coeff of 1st order \nabla = 1/h**2
    e = 10
    L  = np.array(optimal_L_1(n).toarray()) # is an issue if these are not normal arrays
    L_2 = np.array(optimal_L(n).toarray())  # ditto
    #create GEKKO array of variables
    u = m.Array(m.Var,(n**2))
    w = m.Array(m.Var,(n**2))
    #evaluate the variables in arrays
    for i in range(n**2):
        u[i].value = 0.1*ri(0,10) # evenly distibuted points between 0 and 1
        w[i].value = 0.2+(0.01*ri(0,10)) #give small variation to water
    print(u)
    #define u equations  (have to be defined individually rather than vector notation??)                                      
    m.Equations([u[i].dt() == ((u[i]**2)*(w[i])) - (B*u[i]) + dis2*(sum(L_2[i]*u)) for i in range(len(u))])
    
    #define w equations
    m.Equations([0 == A - w[i] - ((u[i]**2)*(w[i])) + dis2*e*(sum(L_2[i]*w))for i in range(len(w))])

    #define time step
    m.time = np.linspace(0,50,100)
    
    #solve DAE
    m.options.IMODE = 5
    m.solve(disp = False)

    print('solved')
    #plot results (copy and pasted from usual solver)
    umax = max([u[i][-1] for i in range(len(u))])
    fig,ax = plt.subplots()
    #for each point
    for x in range(n):
        for y in range(n):
            nnn = u[x+(n*y)][-1]/umax
            if nnn < -1e-2 or nnn>1: ###~~ some small numerical errors in calculation occurred in some instances
                col = 'red'
                print(nnn)
            elif abs(nnn) <= 1e-2:# this covers case of computational numerical errors
                col = 'white'
            
            else:                    
                col = matplotlib.colors.to_hex([1-nnn,1-0.5*nnn,1-nnn])
                    
            ax.plot(x,y, color = col, marker = 's', markersize = 15)
    #disappear surrounding things to look nicer
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_aspect('equal',adjustable='box')
    for item in [fig, ax]:
        item.patch.set_visible(False)
    plt.title('late times results')
    fig.show()
    return u
