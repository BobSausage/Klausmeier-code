import numpy as np
from scipy.integrate import odeint as oi
from scipy.integrate import solve_ivp as si
import matplotlib.pyplot as plt
t = np.linspace(0,20)

#Non-Spatial Klausmeier
#define  z = [ u , w ]
def nskm(z,t, A = 0.5, B = 0.2):
    u = z[0]
    w = z[1]
    dudt = w*(u**2) - B*u
    dwdt = A - w - w*(u**2)
    return [dudt,dwdt]

def ivp_nskm(t,z,A=1.1776,B=0.45):
    u = z[0]
    w = z[1]
    dudt = w*(u**2) - B*u
    dwdt = A - w - w*(u**2)
    return [dudt,dwdt]

#init cond
zo = [1.5,0.8]#[0.8,0.1]

#time same as above

###solve ODE
##z = oi(nskm,zo,t)
##print(z)
###alternatively split into u and w
##u = z[:,0]
##w = z[:,1]
##
###plotting
##########plt.plot(t,z)
##plt.plot(t,u,'g-')
##plt.plot(t,w,'b-')
##plt.xlabel('time')
##plt.ylabel('f(t)')
##plt.legend(['u(t)','w(t)'])
##plt.show()

def roab(amn,amx,bmn,bmx,astep = 0.1, bstep=0.1):
    '''arguments Amin , Amax , Bmin , Bmax, Astep, Bstep,
calculates and shows plot for each value in the range
(if min = max only plots 1 plot in that variable)'''
    for a in range(int(10*amn),int(10*amx)+1):
        for b in range(int(10*bmn),int(10*bmx)+1):
            #solve ODE
            z = oi(nskm,zo,t,(a/10,b/10))
            #define u and w
            u = z[:,0]
            w = z[:,1]

            plt.plot(t,u,'g-')
            plt.plot(t,w,'b-')
            plt.xlabel('time')
            plt.ylabel('f(t)')
            plt.legend(['u(t)','w(t)'])
            plt.title('A = '+str(a/10)+' , B = '+str(b/10))
            plt.show()
                   

def roz(uinit,uend,winit,wend,ustep = 0.1,wstep = 0.05,A=1.1176,B=0.45,pltshow=True):
    '''runs over a range of initial conditions for fixed A and B'''
    outu = [[],[]]
    outw = [[],[]]
    fig, ax = plt.subplots()
    for icu in range(int((1/ustep)*uinit),int((1/ustep)*uend)+1):
        for icw in range(int((1/wstep)*winit),int((1/wstep)*wend)+1):
            #solve ODE
            z = oi(nskm,[(ustep*icu)+uinit,(wstep*icw)+winit],t,(A,B))
            #define u and w
            u = z[:,0]
            w = z[:,1]

            
            
            ax.plot(t,u,'g-')
            ax.plot(t,w,'b-')
            ax.set_xlabel('time')
            ax.set_ylabel('f(t)')
            fig.legend(['u(t)','w(t)'])
            fig.suptitle('u = '+str(round(icu*ustep,3))+' , w = '+str(round(icw*wstep,3)))
            if u[-1] > 1e-2:
                outu[0].append(u[0])
                outw[0].append(w[0])
            else:
                outu[1].append(u[0])
                outw[1].append(w[0])
            if pltshow == True:
                plt.show()
            else:
                fig.clear()
    plt.plot(outu[0],outw[0],'g',marker = 's',markersize = 4,linestyle = 'None')
    plt.plot(outu[1],outw[1],'r',marker = 's',markersize = 4,linestyle = 'None')
    plt.xlabel('u')
    plt.ylabel('w')
    plt.title('The validty of late time plant growth')
    plt.show()



def multi_plot(zlist,a=1.1776, b = 0.45,tmax = 50):
    '''plots different initial conds next to eachother
zlist is a 1 dimensional list of tuples giving i.c.s for u and w'''

    #first define the shape
    if (len(zlist) % 4 == 0) and (len(zlist) > 8):
        #set col = 4,
        col, row = 4 , int(len(zlist)/4)
    elif (len(zlist) % 3 == 0) and (len(zlist) > 6):
        col, row = 3, int(len(zlist)/3)
    elif (len(zlist) % 2 == 0) and (len(zlist) > 3):
        col, row = 2, int(len(zlist)/2)
    else:
        col, row = 1, len(zlist)
        
    #get empty list for axes
    axes = []
    #now define sets of axes
    fig, axes = plt.subplots(row,col)
    #set time span
    tsp = np.linspace(0,tmax,200)
    
    #now solve and plot each plot
    for plot_i in range(col):
        for plot_j in range(row):
            #solve
            semi_sol = si(ivp_nskm,(0,tmax),zlist[plot_i + (plot_j*col)],t_eval=tsp)
            #seperate sols for water and plant
            u = semi_sol.y[0]
            w = semi_sol.y[1]            
            #plot solns
            #try to plot in 2d (this will not work if you want only one row)
            try : 
                axes[plot_j][plot_i].plot(tsp,u,'g')
                axes[plot_j][plot_i].plot(tsp,w,'b')
                axes[plot_j][plot_i].margins(0.02)
            #if plotting 1d then will raise error in which case do this instead
            except:
                axes[plot_j].plot(tsp,u,'g')
                axes[plot_j].plot(tsp,w,'b')
                axes[plot_j].margins(0.02)
    fig.canvas.print_png('Project_plots\multiplot.png')
    plt.show()
            
def basic_fig_plot(z,tmax = 100,a=1.1776,b=0.45,save = False):
    '''z should be (u,w) in form tuple,list or np array'''
    #points it will evaluate at (these are not all points it is evaluated at just the ones you want out)
    #explicitly giving it time points to get out gives a smoother graph
    tsp = np.linspace(0,tmax,1000)
    #solve
    sola = si(ivp_nskm,(0,tmax),z,t_eval = tsp)
    #as z was a list length 2 where u = z[0] and w = z[1] can read solns like this
    u = sola.y[0]
    w = sola.y[1]
    #plot against time points
    plt.plot(tsp,u,'g')
    plt.plot(tsp,w,'b')
    plt.margins(0.01)
    plt.xlabel('time')
    plt.ylabel('density')
    plt.legend(['u(t)','w(t)'])
    plt.rc('axes', labelsize=20)
##    plt.xticks(fontsize = 15)
##    plt.yticks(fontsize = 15)
    
    if save:
        plt.savefig('plots/U_v_W_z' + str(z) + '.png')
    plt.show()
    
