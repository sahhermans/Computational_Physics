# Molecular dynamics simulation Argon - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 31/03/2020

# Datablocking script
import numpy as np
def datablocker(x):
    """
    Returns the standarddeviation of a variable which is determined using
    datablocking. This is done by first splitting the number of variables
    into blocks of equal size. The blocksize has to be chosen in such a way
    that the blocks are no longer correlated. This is achieved when the 
    standarddeviation is more or less constant. 
    
    Input parameters:
        x (float) : an array containing the the variable of which you want
                    to determine the standarddeviation. Note, that the size 
                    of x must be such that it can be split up in equal block-
                    sizes.
                    
    Output:
        sdv (float) : a float that contain the standarddeviation of the input
                      parameter
                        
    
    """

    N=np.size(x)
    blocksize=int(0.04*N)
    
    
    Nb=int(N/blocksize)
    
    a=np.zeros(Nb)
    for b in range(Nb):
        
        a[b]=np.mean(x[b*blocksize:((b+1)*blocksize)])
        
    sdv=np.sqrt(1/(Nb-1)*(np.mean(a**2)-np.mean(a)**2))
    
    return sdv

