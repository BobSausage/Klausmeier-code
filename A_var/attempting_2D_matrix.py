##trying to generate a kronecker matrix.
##simplifies 2D 2nd order diff to 1 matrix removing the for loop

import numpy as np
from scipy import sparse



def getL(n):
        '''gives a sparse matrix L defining discrete 2nd order differentials with Dirichlet b.c.s'''
        L =  -2*np.eye(n)
        L += np.diag(np.ones(n-1),-1)
        L += np.diag(np.ones(n-1),1)
        L[0][n-1] = L[n-1][0] = 1

        #make matrix sparse
        Ls = sparse.csr_matrix(L)

        return Ls
###aim is to create v. large numpy array matrix then squish into sparse matrix


def amat(n):
    mat =  -4*np.eye(n)
    mat += np.diag(np.ones(n-1),-1)
    mat += np.diag(np.ones(n-1),1)
    mat[0][n-1] += 1
    mat[n-1][0] += 1

    return mat

def bmat(n):
    mat = np.eye(n)
    return mat

def omat(n):
    return np.zeros((n,n))



####this only works for n>2 cos realistically never going to need it for n=2 and am lazy
def big_fat_L(n):
    #we want to make columns then append each row together.

    #generate first col of matrices
    #start with diff part
    col_i = amat(n)
    #add identity
    col_i = np.concatenate((col_i,bmat(n)))
    #add n-3 sets of zero
    for q in range(n-3):
        col_i = np.concatenate((col_i,omat(n)))
    #finally finish with another identity
    col_i = np.concatenate((col_i,bmat(n)))

    count = 0
    #now for the middle cols
    #for each col
    for p in range(1,n-1):
        col_p = -1000*np.ones((n,n))## this is just so there is a base it will be removed later
        #for each section of the col
        for q in range(n-2):
            #depending on where col is
            if q ==count:
                #add diagonal
                col_p = np.concatenate((col_p,bmat(n),amat(n),bmat(n)))
            else:
                col_p = np.concatenate((col_p,omat(n)))

        ####now attach this to col_i (this is messy but oh well)
        col_ik = list(col_i)
        for row in range(n**2):
            col_ik[row] = np.concatenate((col_i[row],col_p[row+n]))
        col_i = np.array(col_ik)
##        print(col_i)
        #set ready for next col
        count += 1
##    print(col_i)
    #finally add final column
    col_f = bmat(n)
    for p in range(n-3):
        col_f = np.concatenate((col_f,omat(n)))
    col_f = np.concatenate((col_f,bmat(n),amat(n)))
    
    #finally finally finally attach this column to the matrix
    col_ik = list(col_i)
    for row in range(n**2): 
        col_ik[row] = np.concatenate((col_i[row],col_f[row]))
    col_i = np.array(col_ik)

    return sparse.csr_matrix(col_i)
