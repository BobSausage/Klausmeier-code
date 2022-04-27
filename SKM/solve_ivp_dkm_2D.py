#solveivp
#should behave exactly like dicrete_klaumeier_2D_fn but uses ivp_solve instead to make 

from optimal_L import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import solve_ivp as si
from scipy import sparse
from random import randint as ri
import time



class solve:

    def __init__(self,n,t=11,tmin=0,tmax = 10,A=1.1776,B=0.45,d = 0.25):
        '''n is number of points in each of x and y
t is the number of time poitns,
A&B are the parameneters defined in the eqn'''
        self.n = n
        self.A = A
        self.B = B
        self.t = np.linspace(tmin,tmax,t)
        self.tr = (tmin,tmax)  ##range of t vals (needed for solve_ivp arg)
        self.zsols = np.array([])
        self.z0 = np.array([])##i.c.s empty by default
        self.L = optimal_L(n)
        self.save = 0  # decides whether the plot is saved or not
        self.speed = 0 # 1 = plot without showing to allow quick creation of images
        self.d = d # multiplier of plant dispersal = 1/x**2

    #makes jacobian for certain solvers
    def getjac(self,t,z,e=0):
        '''n is number of points (equivalent to n**2 in main code),
    z is array len(z) = 2n for values of points
    B is parameter in model
    L is the nxn grad**2 operator
    returns 2nx2n sparse matrix which is jacobian evaluated at a point
    e allows for dispersal in water term default is 0'''
        #redefine for ease of reading
        n , B , L = self.n**2 , self.B , self.L
        #define u and w 
        u, w = z[:n] , z[n:]
        #define each quadrant
        fu  = sparse.diags(2*u*w - np.ones(n)*B)
        fuL = fu + self.d*L
        #
        fw  = sparse.diags(u**2)
        #
        gu  = sparse.diags(-2*w*u)
        #
        gw  = sparse.diags((-1)*(np.ones(n) + u**2))
        gwL = gw + (e*L)
        #join pairs of rows (hstack = horizontal stack)
        toptwo = sparse.hstack((fuL,fw))
        bottwo = sparse.hstack((gu,gwL))
        #join for final mat (vstack = vertical stack)
        finp   = sparse.vstack((toptwo,bottwo))
        #make sure is csr as diag sometimes doesn't interact well
        out = sparse.csr_matrix(finp)
        return out
    

    #define the model
    def fn(self,t,z):
        n = self.n
        #first n entries of z = u
        u = np.array(z[:n**2])
        #second n entries of z = w
        w = np.array(z[n**2:])

        #u equation
        du = w*(u**2) -  self.B*u + (self.L*u)*self.d

        #w equation
        dw = self.A*np.ones(n**2) - w - w*(u**2)

        #recombine
        dz = np.concatenate((du,dw))

        return dz

    
    #Solve ODEs
    def solve2d(self, meth = 'LSODA'):
        init_time = time.time()
        #solves ODE
        #######params (func, range of t (2-tuple), points t is evaluated at, method)
        ###########methods are :  nonstiff:{{'RK45', 'RK23' , ' DOP853'}}, stiff:{{'Radau','BDF'}}, general{{LSODA}}
        ##doesn't like being given jac if not necessary
        if meth in ['Radau','BDF']:
            sols = si(self.fn, self.tr,self.z0, t_eval = self.t, method = meth, jac = self.getjac)
        else:
            sols = si(self.fn, self.tr,self.z0, t_eval = self.t, method = meth)
        #get solutions
        self.oggle = sols.t
        self.zsols = sols.y
        #return max value of u & w for all time
        self.umax = max([max(self.zsols[:self.n**2][c]) for c in range(self.n**2)])
        self.wmax = max([max(self.zsols[self.n**2:][c]) for c in range(self.n**2)])

        print('time taken = ' + str(time.time()-init_time))
        return

    
    #shows graph of water and plant densities over time at a given point
    def graph(self,ic):
        '''takes ic tuple (u0,w0)'''
        #split into u solutions and w solutions
        us = self.zsols[ic[0]+(ic[1]*self.n)]
        ws = self.zsols[ic[0]+(ic[1]*self.n)+(self.n**2)]
        #plot
        plt.plot(self.t,us,color = 'g')
        plt.plot(self.t,ws,color = 'b')
        plt.legend(['u(t)','w(t)'])
        plt.show()

    #try and visualise a 2d plot at a given time t
    def visualise(self,t,ms = 4,loc = 'Presentation_plots\giffold'):
        '''t is the ith time interval
ms is the marker size (recommended 4 for 50 points 2 for 150
loc is the location files save to'''
        umax = self.umax
#split plot into constituent parts
        fig,ax = plt.subplots()
        
        #for each point
        for x in range(self.n):
            for y in range(self.n):
                nnn = self.zsols[x+(self.n*y)][t]/umax
                if nnn < -1e-2 or nnn>1: ###~~ some small numerical errors in calculation occurred in some instances
                    col = 'red'
                    print(nnn)
                elif abs(nnn) <= 1e-2:# this covers case of computational numerical errors
                    col = 'white'
            
                else:                    
                    col = matplotlib.colors.to_hex([1-nnn,1-0.5*nnn,1-nnn])
                    
                ax.plot(x,y, color = col, marker = 's', markersize = ms)

        #disappear surrounding things to look nicer
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_aspect('equal',adjustable='box')
        
        for item in [fig, ax]:
            item.patch.set_visible(False)
        if self.save == 1:
            fig.canvas.print_png(loc + '\plotu_t'+ str(t) +'.png')
        if self.speed == 1:
            plt.clf()#clear plot to allow to continue plotting immediately.
            plt.close()
        else:
            plt.title(' t = '  + str(round(self.t[t],2)))
            fig.show()
            


##I.c. reader. Converts a 1D numpy array saved in a txt file to a string

def read_ic(file_name):
    '''reads a file and converts it into a np array
FILE HAS TO CONTAIN A 1-DIM NP ARRAY FORMAT
Must include file extension in file_name'''
    #open file
    file_op   = open(file_name,'r+')
    #read file (put it in a string)
    file_str  = file_op.read()
    #split into lists to remove front and end
    file_lisa = file_str.split('[')
    file_lisb = file_lisa[1].split(']')
    #now we have a string of just numbers and commas
    #split into individual numbers
    file_fin  = file_lisb[0].split(',')
    #make each item a float
    flo_list  = [float(item) for item in file_fin]
    #turn into an np array
    out = np.array(flo_list)
    #sanity check (length should be a square number for i.c.s) so if this is n the ics are compatible
    print(len(out)**0.5)
    #close file
    file_op.close()
    return out























