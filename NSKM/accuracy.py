#added accuracy to try and reduce boundary errors

import numpy as np
from mpmath import fp, mp, mpf, nstr
import matplotlib.pyplot as plt
#mpf = mp float
# nstr is change output digits without changing accuracy of calculation e.g. nstr(a,8) is print 8 dp 

#change precision
##mp.prec = 100 ### sets to 100 bit
#changes decimal places
mp.dps = 1000



#redefine functions same as before but with mpfs
colours = ['red','red','yellow','lime','blue']
#colours : NA  , 1 stbl, C stbl, all stbl, D stable

amin = 0#1.25#1.2
amax = 2#1.28#1.8
bmin = 0#0.62#0.6
bmax = 1#0.64#0.9
pts  = 200#00
save = 1#1

alist = np.linspace(amin,amax,pts)
blist = np.linspace(bmin,bmax,pts)

def getD(A,B):
    a , b = mpf(A), mpf(B)
    return 0.5*((a) - ((a**2 - 4*(b**2))**0.5))

def getC(A,B):
    a , b = mpf(A), mpf(B)
    return 0.5*((a) + ((a**2 - 4*(b**2))**0.5))

for a in alist:
    for b in blist:
        q = 0
        A = mpf(a)
        B = mpf(b)
        if A >= 2*B:
            q = 1
            D = getD(A,B)
            C = getC(A,B)
            evaltc = (B**3 - C**2 - B**2)
            evaldc = (C**2 * B**3) - B**5
            evaltd = (B**3 - D**2 - B**2)
            evaldd = (D**2 * B**3) - B**5
            if (evaltc < 0) and (evaldc > 0):
                q = 2
            if (evaltd < 0) and (evaldd > 0):
                if q == 2:
                    q = 3
                else:
                    q = 4
##            print(A,B,C,D,evaltc,evaldc,evaltd,evaldd)
        plt.plot(A,B,color = colours[q],marker='s',markersize = 1)

plt.margins(0)
plt.xlabel('A', fontsize = 15)
plt.ylabel('B', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
#plt.title('Showing the values of A and B that have either 1,2 or 3 stable points')
if save == 1 :
    res = 720
    #save name encodes data on amin, amax, bmin, bmax and no of points in each dimension
    plt.savefig('plots\plotac_A' + str(amin) + '_' + str(amax) + '_B' + str(bmin) + '_' +str(bmax) + '_' + str(pts)+'.png',dpi=res)
plt.show()


def getpoints(A,B):
    '''takes A and B gives saddle point and stable point'''
    D = getD(A,B)
    C = getC(A,B)
    print(' the saddle point is at (' + str(round(D/B,3)) +','+str(round(C,3))+')')
    print(' the stable point is at (' + str(round(C/B,3)) + ','+str(round(D,3)) + ')')
    
