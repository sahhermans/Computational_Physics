# Molecular dynamics simulation Argon - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 31/03/2020

# Plotfunction
import matplotlib.pyplot as plt
from numpy.random import normal as normal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

def plotfunction(x,y,figurelabel,figurenumber,xlabel,ylabel):
    """
    Makes a plot of given input data
    
    Input parameters:
        x (float)            : array containing all data for the horizontal axis
        y (float)            : array containing all data for the vertical axis
        figurelabel (string) : a label for the data
        figurenumber (int)   : denotes the figurenumber
        xlabel (string)      : provides a label for the horizontal axis
        ylabel (string)      : provides a label for the vertical axis
        
    Output parameters:
        fig (figure)         : a figure of the input data
    """
    
    fig=plt.figure(figurenumber)
    plt.plot(x,y,label= figurelabel)
    plt.xlabel(xlabel, fontsize= 16)
    plt.ylabel(ylabel, fontsize =16)
    #plt.legend(fontsize =13)
    
    return fig

def animation3d(xv,L):
    """
    Makes a 3D animation of the particles.
    
    Input parameters:
        xv (float) : an array containing all the x,y,z compoents of the particles
                     in the simulation
        L (float)  : the length of the box you are going to animate
        
    Output parameters:
        ani  (figure)   : an animation of the dynamics of the particles
    """

    nfr = 4001 # Number of frames
    fps = 10000 # Frame per sec
    xs = xv[:,:,0]
    ys = xv[:,:,1]
    zs = xv[:,:,2]

    fig = plt.figure(5)
    ax = fig.add_subplot(111, projection='3d')
    sct, = ax.plot([], [], [], "o", markersize=6)
    def update(ifrm, xa, ya, za):
        sct.set_data(xa[ifrm], ya[ifrm])
        sct.set_3d_properties(za[ifrm])
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_zlim(0,L)
    ani = animation.FuncAnimation(fig, update, nfr, fargs=(xs,ys,zs), interval=1/fps)
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    #fn = 'plot_3D_scatter_animation_6particles_knownvelpos_v2'
    #ani.save(fn+'.mp4',writer=writer)
    
    return ani

def histogram(x,y,figurelabel,xlabel,ylabel,figurenumber):
    """
    Makes a histogram of the input data. 
    
    Input parameters:
        x (float)            : array containing the data to be binned
        y (int)              : integer giving the amount of bins of the histogram
        figurelabel (string) : a label for the data
        figurenumber (int)   : denotes the figurenumber
        xlabel (string)      : provides a label for the horizontal axis
        ylabel (string)      : provides a label for the vertical axis
        
    Output parameters:
        fig (figure)         : the histogram of the input data
    """
    
    
    fig=plt.figure(figurenumber)
    plt.hist(x, bins = y, label = figurelabel)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    return fig

