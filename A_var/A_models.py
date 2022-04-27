#plotting A simply

import numpy as np
import matplotlib.pyplot as plt

#for model 1
def A(t,Q,R,k,Am):
    return Am+ (k*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R)))

def plot_A(Q=1.3518,R=1,k=1,Am=0.7,ar = (0,8),pltsv=False,pltoff=False):
    b = np.linspace(ar[0],ar[1],100000)
    a = np.array([A(qq,Q,R,k,Am) for qq in b])
    aeff = sum(a)/100000
    if pltoff:
        return aeff
    fig,ax = plt.subplots()
    print('effective A : ',aeff)
    ax.plot(b,a)
    ax.margins(0.0025)
    ax.set_xlabel('t')
    ax.set_ylabel('A(t)')
    ax.set_xbound(0,4)
    ax.set_ybound(0,2.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if pltsv:
        fig.canvas.print_png('plots\\Rainfall_Model_1.png')
    plt.show()
    
def A2(t,Q,R,k1,k2,Am):
    return Am+ (k1*(1 + np.tanh((Q*np.sin(0.5*np.pi*t)) -R))) + (k2*(1 + np.tanh((Q*np.sin((0.5*np.pi*t)+(np.pi))) -R)))

def plot_A2(Q=4,R=3.42,k1=1,k2=0.5,Am=0.7,ar = (0,8),pltsv=False,pltoff=False):
    b = np.linspace(ar[0],ar[1],10000)
    a = np.array([A2(qq,Q,R,k1,k2,Am) for qq in b])
    aeff = sum(a)/10000
    if pltoff:
        return aeff

    fig,ax = plt.subplots()
    print('effective A : ',aeff)
    ax.plot(b,a)
    ax.plot(0,0)
    ax.margins(0.0025)
    ax.set_xlabel('t')
    ax.set_ylabel('A(t)')
    ax.set_xbound(0,4)
    ax.set_ybound(0,2.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if pltsv:
        fig.canvas.print_png('plots\\Rainfall_Model_2.png')
    plt.show()

#####################################################################################
#~ Shell in-/output ~#
######################################################################################

##min_dist = 1000000
##for q in np.linspace(0,10,21):
##	for r in np.linspace(0,10,21):
##		a = plot_A(q,r,1,0.7,pltoff = True)
##		if abs(a-1.1776) < min_dist:
##			Q_out = q
##			R_out = r
##			min_dist = abs(a-1.1776)
##print(Q_out,R_out)
##7.368421052631579 5.263157894736842
## ##RATIO!!!!!!!!! precise value depensd on how curvy the graph is
##
##    0.7174721189591078 
##    0.7288659793814434  v square graph
##    0.7287128712871287
##    0.7529785544082606 v curvy graph

##    maxi = 10000
##    for a in np.linspace(1,1.5,51):
##	qq = plot_A(a,1,1,0.7,pltoff=True)
##	if abs(qq-1.1776) < maxi:
##		maxi = abs(qq-1.1776)
##		out_a = a
##>>> for a in np.linspace(1.35,1.36,501):
##	qq = plot_A(a,1,1,0.7,pltoff=True)
##	if abs(qq-1.1776) < maxi:
##		maxi = abs(qq-1.1776)
##		out_a = a
##
##		
##>>> out_a
##1.3563



##>>> maxi = 10000
##>>> for r in np.linspace(3.41,3.42,101):
##	qq = plot_A2(4,r,1,0.5,pltoff=True)
##	if abs(qq - 1.1776) < maxi:
##		maxi = abs(qq-1.1776)
##		out_r = r
##
##		
##>>> out_r
##3.414
    
    
