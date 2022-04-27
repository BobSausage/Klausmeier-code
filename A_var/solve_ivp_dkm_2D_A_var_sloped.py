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
from statistics import mean


class solve:

    def __init__(self,n,t=11,tmax = 10,Amin=0.7,q=1.3518,r=1,dist = 2,d=100,e=5):
        '''n is number of points in each of x and y
t is the number of time poitns,
A&B are the parameneters defined in the eqn'''
        self.n = n
        self.A = A
        self.B = 0.45
        self.t = np.linspace(0,tmax,t)
        self.tr = (0,tmax)  ##range of t vals (needed for solve_ivp arg)
        self.zsols = np.array([])
        self.z0 = np.array([])##i.c.s empty by default
        self.L_1 = optimal_L_1(n) ##first order laplace
        self.L_2 = optimal_L_2(n) ##second order laplace
        self.save = 0  # decides whether the plot is saved or not
        self.speed = 0 # 1 = plots without showing to allow quick creation of images
        self.dist_1 = 1/(2*dist) # multiplier of plant dispersal = 1/2h_x
        self.dist_2 = 1/(dist**2) # multiplier of plant dispersal = 1/h_x^2
        self.d = d #parameter in model
        self.e = e #parameter in model
        ##allows variable A
        self.al = [] ##list to keep track of values of A at given times
        self.amin = Amin ##params for rainfall function
        self.ka = 1   ##main k for rainfall fn, is k in A1 and k1 in A2
        self.kb = 0.5 ##secondary k only used in A2
        self.q = q   ##param in rainfall fn
        self.r = r    ##rainfall param
        

    #makes jacobian for certain solvers
    def getjac(self,t,z):
        '''n is number of points (equivalent to n**2 in main code),
    z is array len(z) = 2n for values of points
    B is parameter in model
    L is the nxn grad**2 operator
    returns 2nx2n sparse matrix which is jacobian evaluated at a point
    e allows for dispersal in water term default is 0'''
        #redefine for ease of reading
        n , B , d, e   = self.n**2 , self.B , self.d, self.e 
        D1, D2, L1 , L2 = self.dist_1, self.dist_2, self.L_1 , self.L_2
        #define u and w 
        u, w = z[:n] , z[n:]
        #define each quadrant
        fu  = sparse.diags(2*u*w - np.ones(n)*B)
        fuL = fu + D2*L2
        #
        fw  = sparse.diags(u**2)
        #
        gu  = sparse.diags(-2*w*u)
        #
        gw  = sparse.diags((-1)*(np.ones(n) + u**2))
        gwL = gw + (d*D1*L1) + (e*D2*L2)
        #join pairs of rows (hstack = horizontal stack)
        toptwo = sparse.hstack((fuL,fw))
        bottwo = sparse.hstack((gu,gwL))
        #join for final mat (vstack = vertical stack)
        finp   = sparse.vstack((toptwo,bottwo))
        #make sure is csr as diag sometimes doesn't interact well
        out = sparse.csr_matrix(finp)
        return out

    #rainfall fn 1
    def get_A1(self,t):
        Q,R,k,Am = self.q,self.r, self.ka,self.amin
        at = Am+ (k*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R)))
        #average over first period to give an approx a val
        if 0<=t<=4:
            self.al.append(at)
        return at

    #rainfall fn 2
    def get_A2(self,t):
        Q,R,k1,k2,Am = self.q,self.r, self.ka,self.kb,self.amin
        at = Am+ (k1*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R))) + (k2*(1 + np.tanh((Q*np.sin((0.5*np.pi*t)+(np.pi))) -R)))
        #average over first period to give an approx a val
        if 0<t<4:
            self.al.append(at)
        return at

    #define the model
    def fn(self,t,z,rainfall=1):
        n, d, e = self.n, self.d, self.e
        L1, L2, D1, D2 = self.L_1, self.L_2, self.dist_1,self.dist_2

        #doesn't like class fns as args so set here
        if rainfall == 1:
            rainfall = self.get_A1
        else:
            rainfall = self.get_A2
        
        #first n entries of z = u
        u = np.array(z[:n**2])
        #second n entries of z = w
        w = np.array(z[n**2:])
        
        #u equation
        du = w*(u**2) -  self.B*u + (L2*u)*D2
        
        #w equation
        dw = rainfall(t)*np.ones(n**2) - w - w*(u**2) + d*D1*L1*w + e*D2*L2*w
        
        #recombine
        dz = np.concatenate((du,dw))

        return dz

    
    #Solve ODEs
    def solve2d(self, rainfall = 1, meth = 'LSODA'):
        init_time = time.time()
        #solves ODE
        #######params (func, range of t (2-tuple), points t is evaluated at, method)
        ###########methods are :  nonstiff:{{'RK45', 'RK23' , ' DOP853'}}, stiff:{{'Radau','BDF'}}, general{{LSODA}}
        ##doesn't like being given jac if not necessary
        if meth in ['Radau','BDF']:
            sols = si(self.fn, self.tr,self.z0, t_eval = self.t,args = [rainfall], method = meth, max_step = 1e-3, jac = self.getjac)
        else:
            sols = si(self.fn, self.tr,self.z0, t_eval = self.t,args = [rainfall], method = meth, max_step = 1e-3)
        #get solutions
        self.oggle = sols.t
        self.zsols = sols.y
        #return max value of u & w for all time
        self.umax = max([max(self.zsols[:self.n**2][c]) for c in range(self.n**2)])
        self.wmax = max([max(self.zsols[self.n**2:][c]) for c in range(self.n**2)])

        #this is a poor way of approximating a but generally agrees with numerical integration to 3sf
        print('effective A=' + str(mean(self.al))) 
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
                nnn = self.zsols[x+(self.n*y)][t]/self.umax
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
            #makes a trimmed box around graph to get rid of excessive whitespace
            #on some computers this crops too much.
            #can either try changing bounds or just comment out things containing 'bbox' and uncoment last line
            bbox = matplotlib.transforms.Bbox([[0.05,-0.1],[0.9,1]])#([[0,-8],[45,55]]) #second arg works for uni computers for 50x50
            bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(loc + '\plotu_t'+ str(t) +'.png', bbox_inches = bbox)
            #fig.canvas.print_png(loc + '\plotu_t'+ str(t) +'.png')
        if self.speed == 1:
            plt.clf()#clear plot to allow to continue plotting immediately.
            plt.close()
        else:
            plt.title(' t = '  + str(round(self.t[t],2)))
            fig.show()
            
    def custom_L_1(self,xd = 1,yd=1):
        '''Makes L_1 but with specified coeffs to allow different angles
    Can use -ve or +ve xdir, ydir.
    For standard 0,0 is lowest point, n,n is highest, can be inverted by setting xd, yd -ve'''
        n = self.n
        #
        di_vecs = []
        #just off central
        di_vecs.append(xd*np.concatenate((np.tile((np.concatenate((-1*np.ones(n-1),[0]))),n-1),-1*np.ones(n-1))))
        di_vecs.append(xd*np.concatenate((np.tile((np.concatenate((np.ones(n-1),[0]))),n-1),np.ones(n-1))))
        #implement ones in corners of each of the diagonal matrices
        di_vecs.append(xd*np.concatenate((np.tile((np.concatenate(([1],np.zeros(n-1)))),n-1),[1])))
        di_vecs.append(xd*np.concatenate((np.tile((np.concatenate(([-1],np.zeros(n-1)))),n-1),[-1])))
       
        #identity above/belo central diag
        di_vecs.append(yd*-1*np.ones(n**2 - n))
        di_vecs.append(yd*np.ones(n**2 - n))
        #corner identities
        di_vecs.append(yd*np.ones(n))
        di_vecs.append(yd*-1*np.ones(n))
        #
        #get positions of di_vecs
        di_pos = [1,-1,n-1,1-n,n,-n,(n**2)-n,n - (n**2)]
        #
        #make matrix
        mat = sparse.diags(di_vecs,di_pos)
        #0s on diags are saved so to remove change to a csr sparse matrix
        mat = mat.tocsr()
        
        self.L_1 = mat

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























