import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import solve_ivp as si
from imageio import imread
from matplotlib.transforms import Bbox

num = 201 ## edit number of points that special_plot is evaluated at here

### this model is used for the contours and arrows
def model(x, t=0,a = 1.1776,b = 0.45):
    '''non-spatial klausmeier'''
    u,w = x
    du = w*(u**2) - (b*u)
    dw = a - w - w*(u**2)
    return np.array([du,dw])

###the model has arguments for z and t swapped and solves for arrays for background
def ivp_model(t,z,a=1.1776,b=0.45,no_of_evals = num):
    #split into w and u
    u = z[:no_of_evals]
    w = z[no_of_evals:]
    #system
    du = ((u**2)*w) - (b*u)
    dw = a*np.ones(no_of_evals) - w - ((u**2)*w)
    return np.concatenate((du,dw))
    

def special_plot(xrange=(0,4),yrange=(0,2),no_of_evals = num):
    '''very primitive and potentially not very useful visualisation of where points end up plotted with the contours of the system'''

    outu , outw = [[],[]] , [[],[]]
    no_of_wvals = int((no_of_evals + 1)/2)
    
    ####first we want to plot ~5,000 points and where they are evaluated
    for wi in np.linspace(yrange[0],yrange[1],no_of_wvals):
        #create a row of initial conds. 
        zi = np.concatenate((np.linspace(xrange[0],xrange[1],no_of_evals),wi*np.ones(no_of_evals)))
        #solve in ode solver
        zsolved = si(ivp_model,(0,100),zi)
        #get information we care about
        usols = zsolved.y[0:no_of_evals] ##removes the w coords

        #now determine if each point is +ve or -ve
        for ui in range(no_of_evals):
            #if the plant population at late times is significant then add to the first list
            if usols[ui][-1] > 1e-3:
                outu[0].append(zi[ui])
                outw[0].append(wi)
            else:
                outu[1].append(zi[ui])
                outw[1].append(wi)
    fig, ax = plt.subplots()            
    #now all vals solved and sorted plot them
    ax.plot(outu[0],outw[0],color = matplotlib.colors.to_hex([0.76,1,0.76]),marker = 's',markersize = 2.5,linestyle = 'None')
    ax.plot(outu[1],outw[1],color = matplotlib.colors.to_hex([0.8,0.8,1]),marker = 's',markersize = 2.5,linestyle = 'None')
    ax.set_xlabel('u')
    ax.set_ylabel('w')
    ###gets rid of space between axis and plot
    plt.margins(0)
    ##now cut out the actual plot fromt he axes
    #bbox = Bbox([[xrange[0],yrange[0]],[xrange[1],yrange[1]]])
    bbox = Bbox([[0,0],[1,1]])
    #
    bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())
    #save fig
    fig.savefig('wells.png', bbox_inches = bbox)
    #show progress and clear plot
    plt.show()

    #opens saved figure
    im = imread('wells.png')
    #scales figure to size of plot
    aspect = im.shape[0]/im.shape[1] * (xrange[1]-xrange[0])/(yrange[1]-yrange[0])
    #sets figure to background of plot
    plt.imshow(im, extent = [xrange[0],xrange[1],yrange[0],yrange[1]], aspect = aspect)
    ##Now plot the contours and arrows over top
    x_ = np.linspace(xrange[0], xrange[1], 25)                                                             
    y_ = np.linspace(yrange[0], yrange[1], 25)                                                             

    grid = np.meshgrid(x_, y_)

    dfmat = np.zeros((25,25, 2))
    for nx in range(25):
        for ny in range(25):
            df = model([grid[0][nx,ny], grid[1][nx,ny]])
            dfmat[nx, ny, 0] = 0.00000005*df[0]
            dfmat[nx, ny, 1] = 0.00000005*df[1]

    ###added code below
    #for contours
    x1 = np.linspace(xrange[0], xrange[1], 200)                                                             
    y1 = np.linspace(yrange[0], yrange[1], 200)
    grid2 = np.meshgrid(x1, y1)
    dfmat2 = np.zeros((200,200, 2))
    for nx in range(200):
        for ny in range(200):
            df = model([grid2[0][nx,ny], grid2[1][nx,ny]])
            dfmat2[nx, ny, 0] = df[0]
            dfmat2[nx, ny, 1] = df[1]
    
    #plots arrows we want same no of arrows different no of points in line
    for aaaa in range(1,12):
        for bbbb in range(1,12):
            aaa, bbb = aaaa*2, bbbb*2
            plt.arrow(grid[0][aaa][bbb],grid[1][aaa][bbb],dfmat[:,:,0][aaa][bbb],dfmat[:,:,1][aaa][bbb],head_width=0.025, color = 'black',length_includes_head = True)
    #plots u dot contour
    plt.contour(grid2[0], grid2[1], dfmat2[:, :, 0], [0], colors = 'r')
    #plots w dot contour
    plt.contour(grid2[0], grid2[1], dfmat2[:, :, 1], [0], colors = 'g')
    plt.xlabel('u', fontsize = 15)
    plt.ylabel('w', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.savefig('Project_plots/phaseplane_smooth_A_1_1776__B_0_45.png')
    
    plt.show()

def from_pot(xrange=(0,4),yrange=(0,2),no_of_evals = num):
    
    #opens saved figure
    im = imread('wells.png')
    #scales figure to size of plot
    aspect = im.shape[0]/im.shape[1] * (xrange[1]-xrange[0])/(yrange[1]-yrange[0])
    #sets figure to background of plot
    plt.imshow(im, extent = [xrange[0],xrange[1],yrange[0],yrange[1]], aspect = aspect)
    ##Now plot the contours and arrows over top
    x_ = np.linspace(xrange[0], xrange[1], 25)                                                             
    y_ = np.linspace(yrange[0], yrange[1], 25)                                                             

    grid = np.meshgrid(x_, y_)

    dfmat = np.zeros((25,25, 2))
    for nx in range(25):
        for ny in range(25):
            df = model([grid[0][nx,ny], grid[1][nx,ny]])
            dfmat[nx, ny, 0] = 0.0000005*df[0]
            dfmat[nx, ny, 1] = 0.0000005*df[1]

    ###added code below
    #for contours
    x1 = np.linspace(xrange[0], xrange[1], 200)                                                             
    y1 = np.linspace(yrange[0], yrange[1], 200)
    grid2 = np.meshgrid(x1, y1)
    dfmat2 = np.zeros((200,200, 2))
    for nx in range(200):
        for ny in range(200):
            df = model([grid2[0][nx,ny], grid2[1][nx,ny]])
            dfmat2[nx, ny, 0] = df[0]
            dfmat2[nx, ny, 1] = df[1]
    
    #plots arrows we want same no of arrows different no of points in line
    for aaaa in range(1,12):
        for bbbb in range(1,12):
            aaa, bbb = aaaa*2, bbbb*2
            plt.arrow(grid[0][aaa][bbb],grid[1][aaa][bbb],dfmat[:,:,0][aaa][bbb],dfmat[:,:,1][aaa][bbb],head_width=0.025, color = 'black',length_includes_head = True)
    
    #plots u dot contour
    plt.contour(grid2[0], grid2[1], dfmat2[:, :, 0], [0], colors = 'r')
    #plots w dot contour
    plt.contour(grid2[0], grid2[1], dfmat2[:, :, 1], [0], colors = 'g')
    plt.xlabel('u', fontsize = 15)
    plt.ylabel('w', fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.margins(0)
##    plt.savefig('Project_plots/phaseplane_smooth_A_1_1776__B_0_45.png')
    
    plt.show()


    
if __name__ == "__main__":
    special_plot()
    #from_pot()
