# Molecular dynamics simulation Argon - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 25/03/2020

# Initialization script

import math
import numpy as np
import statistics as stats

def fcc(N,D,a,ncopies,timesteps):
    """ 
    Returns an ffc lattice configuration. Function initializes a primitive unit 
    cell of an FCC lattice and copies this over all dimensions.
    
    Input parameters:
        N (int): gives the amount of particles in the system
        
        D (int): gives the amount of dimensions of the system
        
        a (float): gives the lattice constant of interest
        
        ncopies (int): gives the amount of copies of the unit cell 
                       of the FCC lattice in each dimension.
        
        timesteps (int): gives the total number of timesteps
        
    Returned: 
        xv (float): array containing all x-, y- and z-components
        of all particles in an ffc configuration
    """
    xv = np.zeros((timesteps + 1, N, D), dtype = float)
    punitcell = np.array([[0,0,0],[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0]]) # primitive unit cell of FCC lattice
    particle = 0
    for g in range(len(punitcell)):
        for s in range(ncopies):
            for d in range(ncopies):
                for f in range(ncopies):
                    xv[0,particle,:] = a/4 + punitcell[g] + a * np.array([s,d,f])
                    particle = particle + 1
                    
    return xv

def randompos(L,N,D,timesteps):
    """
    Returns a random distribution of N particles in a D-dimensional box with
    length L. 
    
    Input parameters:
        N (int): gives the amount of particles in the system
        
        D (int): gives the amount of dimensions of the system
        
        L (float): gives the length of the simulation box into each dimension
             
        timesteps (int): gives the total number of timesteps
        
    Returned: 
        xv (float): array containing all x-, y- and z-components
        of all particles in the system
    
    """
    xv = np.zeros((timesteps + 1, N, D), dtype = float)
    
    xv[0,:,:] = L * np.random.rand(N,D)
    
    return xv

def checkpos(timesteps):
    """
    Initializes 6 particles in a known 3-dimensional pattern to enable checking 
    the functioning of the simulation.
    
    !Warning! Should use this function together with the 'checkvel' function
    to initialize the velocities.
   
    Input parameters:
        N (int): there are 6 particles in this system
        
        D (int): the system is three dimensional
        
        timesteps (int): gives the total number of timesteps
        
    Returned: 
        xv (float): array containing all x-, y- and z-components
        of all 6 particles in the system
    """
    N = 6
    D = 3
    
    xv = np.zeros((timesteps + 1, N, D), dtype = float)
    
    xv[0,:,:]=np.array([[1.5,2.5,2.5],[3.5,2.5,2.5],[2.5,1.5,2.5],[2.5,3.5,2.5],[2.5,2.5,1.5],[2.5,2.5,3.5]])

    return xv

def maxwellboltzmann(T,N,D,timesteps):
    """
    Initializes the velocities of N particles according to the 
    Maxwell-Boltzmann distribution. 
    
    For-loop at the end makes sure that the mean velocity of the system 
    is zero.
        
    Input parameters:
        T (float): gives the initialy set temperature in dimensionless
        units
        
        N (int): gives the amount of particles in the system
        
        D (int): gives the amount of dimensions of the system
        
        timesteps (int): gives the total number of timesteps
        
    Returned: 
        vv (float): array containing all x-, y- and z-velocity components
        of all particles in the system
    
    """
    vv = np.zeros((timesteps + 1, N, D), dtype = float)
    vv[0,:,:] = np.random.normal(0,math.sqrt(T),(N,D))

    for i in range(D):
        vv[0,:,i] = vv[0,:,i] - stats.mean(vv[0,:,i])

    return vv

def checkvel(timesteps):
    """   
    Initializes the velocities of the 6 particles from the 'checkpos' function
    in such a way that they all move towards each other.
    
    For-loop at the end makes sure that the mean velocity of the system 
    is zero.
   
    !Warning! Should use this function together with the 'checkpos' function
    to initialize the positions.
    
    Input parameters:
        N (int): there are 6 particles in this system
        
        D (int): the system is three dimensional
        
        timesteps (int): gives the total number of timesteps
        
    Returned: 
        vv (float): array containing all x-, y- and z-velocity components
        of all 6 particles in the system
    
    """
    N = 6
    D = 3
    
    vv = np.zeros((timesteps + 1, N, D), dtype = float)
    
    vv[0,:,:]=np.array([[2,0,0],[-2,0,0],[0,2,0],[0,-2,0],[0,0,2],[0,0,-2]]) # 6 particles moving towards each other

    for i in range(D):
        vv[0,:,i] = vv[0,:,i] - stats.mean(vv[0,:,i])

    return vv