#defining the jacobian function

import numpy as np
from scipy import sparse

def getjac(n,z,B,L,e=0):
    '''n is number of points (equivalent to n**2 in main code),
z is array len(z) = 2n for values of points
B is parameter in model
L is the nxn grad**2 operator
returns 2nx2n sparse matrix which is jacobian evaluated at a point
e allows for dispersal in water term default is 0'''
    #define u and w 
    u, w = z[:n] , z[n:]
    #define each quadrant
    fu  = (2*u*w - np.ones(n)*B)*np.eye(n)
    fuL = fu + L
    fw  = (u**2)*np.eye(n)
    #
    gu  = (-2*w*u)*np.eye(n)
    gw  = (-1)*(np.ones(n) + u**2)*np.eye(n)
    gwL = gw + (e*L)
    #had to keep matrices as np array to allow concatenation
    #join pairs of rows
    toptwo = np.concatenate((fuL,fw),axis=1)
    bottwo = np.concatenate((gu,gwL),axis=1)
    #join for final mat
    finp   = np.concatenate((toptwo,bottwo))
    #convert to sparse
    out = sparse.csr_matrix(finp)
    return out

###new optimised jac getter (around 100* faster on 50x50)
def getjac2(n,z,B,L,e=0):
    '''n is number of points (equivalent to n**2 in main code),
z is array len(z) = 2n for values of points
B is parameter in model
L is the nxn grad**2 operator
returns 2nx2n sparse matrix which is jacobian evaluated at a point
e allows for dispersal in water term default is 0'''
    #define u and w 
    u, w = z[:n] , z[n:]
    #define each quadrant
    fu  = sparse.diags(2*u*w - np.ones(n)*B)
    fuL = fu + L
    fw  = sparse.diags(u**2)
    #
    gu  = sparse.diags(-2*w*u)
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

if __name__ == '__main__':
    from optimal_L import *

    
