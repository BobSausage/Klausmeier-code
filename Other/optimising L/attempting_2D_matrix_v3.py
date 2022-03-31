##trying to generate a kronecker matrix.
##simplifies 2D 2nd order diff to 1 matrix removing the for loop

import numpy as np
from scipy import sparse



def getL(n):
        '''gives a sparse matrix L defining discrete 2nd order differentials with periodic b.c.s'''
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
    #we want to make columns then append to make matrix.
    #first get a sample of each of the 3 types of matrix
    am, bm, om = amat(n),bmat(n), omat(n)
    #generate first col of matrices
    #start with diff part
    col_i = am
    #add identity
    col_i = np.concatenate((col_i,bm))
    #add n-3 sets of zero
    for q in range(n-3):
        col_i = np.concatenate((col_i,om))
    #finally finish with another identity
    col_i = np.concatenate((col_i,bm))

    count = 0
    #now for the middle cols
    #for each col
    for p in range(1,n-1):
        #if this is first one
        if count == 0:
            col_p = np.concatenate((bm,am,bm))
        else:
            col_p = om   
        
        #for each section of the col
        for q in range(1,n-2):
            #depending on where col is
            if q ==count:
                #add diagonal
                col_p = np.concatenate((col_p,bm,am,bm))
            else:
                col_p = np.concatenate((col_p,om))
##        print(col_p)
        ####now attach this to col_i 
        col_i = np.concatenate((col_i,col_p),axis = 1)

        #set ready for next col
        count += 1
##    print(col_i)
    #finally add final column
    col_f = bm
    for p in range(n-3):
        col_f = np.concatenate((col_f,om))
    col_f = np.concatenate((col_f,bm,am))
    #add col-f
    col_i = np.concatenate((col_i,col_f),axis= 1)

    #print(col_i)
    return sparse.csr_matrix(col_i)


#this works for all n>1 and is much faster
def new_L(n):
    #first set diagonals we want diag and either side for in 1 dim
    #central
    mat = -4*np.eye(n**2)
    #off centrals
    #require n-1 1s then a 0 (tile allows this to be repeated n times) (this has to have n**2 - 1 components but some are 0)
    off_pat = np.concatenate((np.tile((np.concatenate((np.ones(n-1),[0]))),n-1),np.ones(n-1)))
    mat += np.diag(off_pat,-1)
    mat += np.diag(off_pat,1)
    #add identity on extra diags for other direction
    mat += np.diag(np.ones((n**2)-n),n)
    mat += np.diag(np.ones((n**2)-n),-n)
    #now add corner identities (for periodic b.c.s in one direction)
    mat += np.diag(np.ones(n),(n**2)-n)
    mat += np.diag(np.ones(n),n - (n**2))
    #finally implement the periodic b.c.s in other direction
    off_pat = np.concatenate((np.tile((np.concatenate(([1],np.zeros(n-1)))),n-1),[1]))
    mat += np.diag(off_pat,n-1)
    mat += np.diag(off_pat, 1-n)
    #make sparse
    return sparse.csr_matrix(mat)

##trying to optimise big_fat_L for larger matrices...
##sparses after concat a row
def biggie2(n):
    #we want to make columns then append to make matrix.
    #first get a sample of each of the 3 types of matrix
    am, bm, om = amat(n),bmat(n), omat(n)
    #generate first col of matrices
    #start with diff part
    col_i = am
    #add identity
    col_i = np.concatenate((col_i,bm))
    #add n-3 sets of zero
    for q in range(n-3):
        col_i = np.concatenate((col_i,om))
    #finally finish with another identity
    col_i = np.concatenate((col_i,bm))
    col_i = sparse.csr_matrix(col_i)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    count = 0
    #now for the middle cols
    #for each col
    for p in range(1,n-1):
        #if this is first one
        if count == 0:
            col_p = np.concatenate((bm,am,bm))
        else:
            col_p = om   
        
        #for each section of the col
        for q in range(1,n-2):
            #depending on where col is
            if q ==count:
                #add diagonal
                col_p = np.concatenate((col_p,bm,am,bm))
            else:
                col_p = np.concatenate((col_p,om))
##        print(col_p)
        ####now attach this to col_i

        col_p = sparse.csr_matrix(col_p)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
        col_i = sparse.hstack((col_i,col_p))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #set ready for next col
        count += 1
##    print(col_i)
    #finally add final column
    col_f = bm
    for p in range(n-3):
        col_f = np.concatenate((col_f,om))
    col_f = np.concatenate((col_f,bm,am))
    #add col-f
    col_f = sparse.csr_matrix(col_f)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    col_i = sparse.hstack((col_i,col_f))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #print(col_i)
    return col_i
        
##sparses after concat a row
def biggie3(n):
    #we want to make columns then append to make matrix.
    #first get a sample of each of the 3 types of matrix
    am, bm, om = sparse.csr_matrix(amat(n)),sparse.csr_matrix(bmat(n)), sparse.csr_matrix(omat(n))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #generate first col of matrices
    #start with diff part
    col_i = am
    #add identity
    col_i = sparse.vstack((col_i,bm))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #add n-3 sets of zero
    for q in range(n-3):
        col_i = sparse.vstack((col_i,om))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #finally finish with another identity
    col_i = sparse.vstack((col_i,bm))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #col_i = sparse.csr_matrix(col_i)

    count = 0
    #now for the middle cols
    #for each col
    for p in range(1,n-1):
        #if this is first one
        if count == 0:
            col_p = sparse.vstack((bm,am,bm))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else:
            col_p = om   
        
        #for each section of the col
        for q in range(1,n-2):
            #depending on where col is
            if q ==count:
                #add diagonal
                col_p = sparse.vstack((col_p,bm,am,bm))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else:
                col_p = sparse.vstack((col_p,om))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        print(col_p)
        ####now attach this to col_i
                
        col_i = sparse.hstack((col_i,col_p))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #set ready for next col
        count += 1
##    print(col_i)
    #finally add final column
    col_f = bm
    for p in range(n-3):
        col_f = sparse.vstack((col_f,om))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    col_f = sparse.vstack((col_f,bm,am))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #add col-f
    col_i = sparse.hstack((col_i,col_f))#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #print(col_i)
    return col_i

####optimised split function
def final_L(n):
    #splitting as new_l is more efficient for small n
    if n < 55:
        return new_L(n)
    else:
        return biggie3(n)




if __name__ == '__main__':
        a = omat(5)
        b = sparse.csr_matrix(a)
        c = amat(5)
        d = sparse.csr_matrix(c)
        e = sparse.vstack((b,d))
        print(e)
        
