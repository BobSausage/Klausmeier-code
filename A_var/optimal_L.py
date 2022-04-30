#neat version of attempting_2D_matrix
#contains 2nd order laplacian approx imported to solve_ivp in solve class

from scipy import sparse
import numpy as np

def optimal_L_2(n):
    #put all vectors along diagonals in a list
    di_vecs = []
    #main vector
    di_vecs.append(-4*np.ones(n**2))
    #just off central
    di_vecs.append(np.concatenate((np.tile((np.concatenate((np.ones(n-1),[0]))),n-1),np.ones(n-1))))
    di_vecs.append(np.concatenate((np.tile((np.concatenate((np.ones(n-1),[0]))),n-1),np.ones(n-1))))
    #identity above/belo central diag
    di_vecs.append(np.ones(n**2 - n))
    di_vecs.append(np.ones(n**2 - n))
    #corner identities
    di_vecs.append(np.ones(n))
    di_vecs.append(np.ones(n))
    #implement ones in corners of each of the diagonal matrices
    di_vecs.append(np.concatenate((np.tile((np.concatenate(([1],np.zeros(n-1)))),n-1),[1])))
    di_vecs.append(np.concatenate((np.tile((np.concatenate(([1],np.zeros(n-1)))),n-1),[1])))
    #
    #get positions of di_vecs
    di_pos = [0,1,-1,n,-n,(n**2)-n,n - (n**2),n-1,1-n]
    #
    #make matrix
    mat = sparse.diags(di_vecs,di_pos)
    #0s on diags are saved so to remove change to a csr sparse matrix
    mat = mat.tocsr()
    
    return mat

def optimal_L_1(n):
    ###it is important to notice L1 is NOT symmetric so careful with ordering
    #put all vectors along diagonals in a list
    di_vecs = []
    #just off central
    di_vecs.append(np.concatenate((np.tile((np.concatenate((-1*np.ones(n-1),[0]))),n-1),-1*np.ones(n-1))))
    di_vecs.append(np.concatenate((np.tile((np.concatenate((np.ones(n-1),[0]))),n-1),np.ones(n-1))))
    #identity above/belo central diag
    di_vecs.append(-1*np.ones(n**2 - n))
    di_vecs.append(np.ones(n**2 - n))
    #corner identities
    di_vecs.append(np.ones(n))
    di_vecs.append(-1*np.ones(n))
    #implement ones in corners of each of the diagonal matrices
    di_vecs.append(np.concatenate((np.tile((np.concatenate(([1],np.zeros(n-1)))),n-1),[1])))
    di_vecs.append(np.concatenate((np.tile((np.concatenate(([-1],np.zeros(n-1)))),n-1),[-1])))
    #
    #get positions of di_vecs
    di_pos = [1,-1,n,-n,(n**2)-n,n - (n**2),n-1,1-n]
    #
    #make matrix
    mat = sparse.diags(di_vecs,di_pos)
    #0s on diags are saved so to remove change to a csr sparse matrix
    mat = mat.tocsr()
    
    return mat
