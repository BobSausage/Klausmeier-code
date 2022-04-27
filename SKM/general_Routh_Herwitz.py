#RH criterion for eps and rho_k
#heavily based on "accuracy" which evaluates RH crit for A B vals

from mpmath import fp, mp, mpf, nstr
import matplotlib.pyplot as plt
import numpy as np

mp.dps = 1000

def getG(eps,rho,A=1.1776,B=0.45):
    ep,ro,a,b = mpf(eps),mpf(rho),mpf(A),mpf(B)
    c  =  0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    g  = 1 + ((c**2)/(b**2)) - (ep*b) + (ep*ro)
    if ep == 0:
        return 0
    return g/ep

def getG_vara(A,rho,eps=1,B=0.45):
    ep,ro,a,b = mpf(eps),mpf(rho),mpf(A),mpf(B)
    if a < 2*b:
        return 0
    c  =  0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    g  = 1 + ((c**2)/(b**2)) - (ep*b) + (ep*ro)
    if ep == 0:
        return 0
    return g/ep
def getH(eps,rho,A=1.1776,B=0.45):
    ep,ro,a,b = mpf(eps),mpf(rho),mpf(A),mpf(B)
    c  = 0.5*(a + (a**2 - (4*b**2))**0.5)
    h  = ((ro - b)*(1 + ((c**2)/(b**2)))) + (2*(c**2)/b)
    if ep == 0:
        return 0
    return h/ep

col = ['black','red','yellow','lime','blue']

def RH_gen(vara_r=(0,1),varb_r=(-1,1),fn_a=getG,fn_b=getH,A=1.1776,B=0.45,points=100):
    '''takes tuple for variable 1 & 2 ranges returns plot of areas they are feasible for A and B'''
    # get points
    e_range = np.linspace(vara_r[0],vara_r[1],points)
    r_range = np.linspace(varb_r[0],varb_r[1],points)
    # run over all points
    for a in e_range:
        for b in r_range:
            #set colour to black
            co = 0
            #evaluate G
            gg = fn_a(a,b,A,B)
            #check RH satisfied
            if gg>0:
                co += 1
            #evaluate H
            hh = fn_b(a,b,A,B)
            #check RH satisfied
            if hh>0:
                co+= 0.5
            #plot point
            plt.plot(b,a,marker = 's',markersize=3.1,color = col[int(2*co)])
    plt.margins(0)
    plt.xlabel('rho')
    plt.ylabel('epsilon')
    plt.show()



#case 2
def get_lambda(rho,e,A,B):
    eh,ro,a,b = mpf(e),mpf(rho),mpf(A),mpf(B)
    c =  0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    num = (ro-b)*(1 + ((c/b)**2) + (eh*ro)) + (2*(c**2)/b)
    den = 1 + ((c/b)**2) + ((eh)*ro)
    return num/den

def get_lam2(rho,e,A,B):
    eh,ro,a,b = mpf(e),mpf(rho),mpf(A),mpf(B)
    c =  0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    num = -2*b*(c**2)
    


def lambda_k(vara_r=(0,1),varb_r=(-1,1),fn_b=get_lambda,A=1.1776,B=0.45,points=100):
    #get points
    a_range = np.linspace(vara_r[0],vara_r[1],points)
    b_range = np.linspace(varb_r[0],varb_r[1],points)
    for a in a_range:
        for b in b_range:
            q = fn_b(a,b,A,B)
            if q > 0:
                col = 'g'
            elif q < 0:
                col = 'r'
            else:
                col = 'y'
            plt.plot(a,b,marker='s',markersize = 3.1, color = col)
    plt.margins(0)
    plt.xlabel('rho')
    plt.ylabel('e')
    plt.show()


def lam_fin(a,b,e,x):
    a,b,e = mpf(a),mpf(b),mpf(e)
    c = 0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    num = 2*(c**2)*b
    den = (b**2) + (c**2) + ((b**2)*e*x)
    out = b - x - (num/den)
    return out

def find_max_rho(a,b,e,x):
    a,b,e = mpf(a),mpf(b),mpf(e)
    c = 0.5*(a + ((a**2)-(4*(b**2)))**0.5)
    num = 2*(c**2)*(b**3)*e
    den = ((b**2) + (c**2) + ((b**2)*e*x))**2
    return num/den

##describes where maxima is found
def rho_v_e(a,b,emax=50,points=201,save = False):
    e_range = np.linspace(0,emax,points)
    rho_ran = np.linspace(0,0.5,points)
    for e in e_range:
        for r in rho_ran:
            ev = find_max_rho(a,b,e,r)
            if ev > 1:
                col = 'r'
            else:
                col = 'yellow'
            plt.plot(r,e,color = col,marker='s',markersize = 3.25)
    plt.margins(0)
    plt.xlabel('rho')
    plt.ylabel('e')
    plt.rc('axes', labelsize=20) 
    if save:
        plt.savefig('Project_plots\\rho_v_e_maxima_a'+str(a)+'__b_'+str(b)+'.png')
    plt.show()
#describes where lambda >0
def rho_v_e_prop(a,b,emax=50,points=201,save = False):
    e_range = np.linspace(0,emax,points)
    rho_ran = np.linspace(0,0.5,points)
    for e in e_range:
        for r in rho_ran:
            ev = lam_fin(a,b,e,r)
            if ev > 0:
                col = 'r'
            else:
                col = 'yellow'
            plt.plot(r,e,color = col,marker='s',markersize = 3.25)
    plt.margins(0)
    plt.xlabel('rho')
    plt.ylabel('e')
    plt.rc('axes', labelsize=20) 
    if save:
        plt.savefig('Project_plots\\rho_v_e_lambda_a'+str(a)+'__b_'+str(b)+'.png')
    plt.show()
