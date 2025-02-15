# Molecular dynamics simulation Argon - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 31/03/2020

# Simulation script

import numpy as np
import math

def force(N,D,L,xv,time,r,Ur,dUdr,F):
    """ 
    Returns the total force on each particle at the time 'time'. To do this 
    the distance between each particle and all others is calculated, as
    well as the derivative of the potential energy for each particle.
    For simplicity, the total potential energy of the system at the time 
    time is also calculated. 
    
    Input parameters:
        N (int): gives the amount of particles in the system
        
        D (int): gives the number of dimensions of the system
       
        L (float): gives the length of the simulation box into each dimension
        
        xv (float): array containing all x-, y- and z-components
        of all particles in the system
        
        time (int): gives the time of interest
        
    Output:
        dv (float): array containing the difference vectors between the
                    particle of interest and all other particles
        
        r (float): array containing all absolute distances between particles
                   r[time,i,j]: absolute distance between particle i and 
                   particle j at time 'time'
        
        Ur (float): array that gives the total potential energy in the system 
        
        dUdr (float): array that gives the derivative of the LJ-potential 
                      for a particle at time 'time'
        
        F (float): returned array containing the force on each particle at time 
                   time.
    """
    dv = np.zeros((N,D), float)
    for i in range(N):
        dv = ((np.tile(xv[time,i,:],(N,1)) - xv[time,:,:]) + L/2) % L - L/2
        r[time,i,:] = np.sqrt(np.sum(dv**2,axis = 1))
        r[time,i,i] = 1
        Ur[time,:] = Ur[time,:] + np.sum((4*((r[time,i,:])**(-12) - (r[time,i,:])**(-6))))/2
        dUdr[time,i,:] = 4*(-12*r[time,i,:]**(-13) + 6*r[time,i,:]**(-7)) 
        dUdr[time,i,i] = 0
                
        F[time,i,:] = np.dot(-dUdr[time,i,:]/r[time,i,:],dv)

    return F[time,:,:]
    
def simulation(N,D,T,L,m,xv,vv,timesteps,equitime,h,nbins):
    """ 
    Performs the molecular dynamics simulation, making use of the force function. 
    
    The function firstly initializes all variables of interest. Subsequently it 
    runs a loop over all timesteps, evaluating the positions and velocities of 
    all particles in the system at each timestep using the Verlet algorithm and
    the force function. 
    
    Simultaneously it calculates different variables of interest, like the virial 
    (total average and for each timestep), the kinetic energy, the total energy 
    and the force. 
    
    Untill the equilibration time equitime, the kinetic energy of the system is 
    rescaled every 10 timesteps by comparing the kinetic energy of the system to 
    the theoretical kinetic energy of the system, determined by the equipartition 
    theory.
    
    After exiting the loop over all timesteps, the function calculates multiple 
    variables of interest, like the temperature, the heat capacity, the 
    compressibility and the pair correlation function.
    
    Input parameters:
        N (int): gives the amount of particles in the system
        
        D (int): gives the amount of dimensions of the system
        
        T (float): gives the initialy set temperature in dimensionless
        units
        
        L (float): gives the length of the simulation box into each dimension
        
        m (float): mass of the simulated particles of interest
        
        xv (float): array containing all x-, y- and z-components
        of all particles in the system
        
        vv (float): array containing all x-, y- and z-velocity components
        of all 6 particles in the system
               
        timesteps (int): gives the total number of timesteps
        
        equitime (int): gives the number of timesteps for the system
                        untill which equilibrates by rescaling of 
                        the kinetic energy
        
        h (float): gives the length of a timestep in unitless time
        
        nbins (int): gives the amount of bins for the pair correlation
                     function
        
    Returned:
    
        Tave (float): the average temperature of the system after
                      equilibration
        
        Taveautocor (float): array containing the temperature of the 
                             system at each timestep, to be used in
                             the autocorrelation function 
        
        Uave (float): the average potential energy per particle in the 
                      system after equilibration
        
        Cave (float): the average heat capacity of the system after
                      equilibration
        
        compressibility (float): the compressibility of the system
        
        compressibilityautocor (float): array containing the compressibility
                                        of the system at each timestep, to 
                                        be used in the autocorrelation 
                                        function
        
        r (float): array containing all absolute distances between particles
                   r[time,i,j]: absolute distance between particle i and 
                   particle j at time 'time'
        
        Ur (float): array that gives the total potential energy in the system
                    at each timestep
        
        dUdr (float): array that gives the derivative of the LJ-potential 
                      for a particle at time 'time'
        
        F (float): array containing the force on each particle at time 
                   time 
    
        xv (float): array containing all x-, y- and z-components
        of all particles in the system
        
        vv (float): array containing all x-, y- and z-velocity components
        of all 6 particles in the system
    
        Ekin (float): array containg the kinetic energy in the system
                      at each timestep
        
        Energytotal (float): array containing the total energy in the system
                             (kinetic + potential) at each timestep
        
        rhist (float): array containing the coordinates of the centers of the 
                       bins of the histogram used to calculte the pair 
                       correlation function
        
        hist (float): array containing the amount of histogram counts in each
                      bin for the interatomic distances between particles at
                      each timestep
        
        gr (float): the time-averaged pair correlation function.
        
        grauto (float): array containing the evaluation of the pair correlation
                        function at each timestep, to be used in the 
                        autocorrelation function.
    """

    # Initialize particle simulation for N particles in D dimensions
    #The extension autocor means that this variable is later used for the autocorrelation function
    r = np.zeros((timesteps+1,N,N), float)
    Ur = np.zeros((timesteps+1,1), float)
    dUdr = np.zeros((timesteps+1,N,N), float)
    Ekin = np.zeros((timesteps+1,1), float)
    F = np.zeros((timesteps+1,N,D), float)
    Energytotal = np.zeros((timesteps+1,1), float)
    virial = 0
    virialautocor = np.zeros([timesteps-equitime-1])
    lam = np.ones(timesteps) * 2
    hist = np.zeros(nbins)
    
    Ekinequi = (3/2)*(N-1)*T # kinetic energy given by equipartition theory
    
    F[0,:,:] = force(N,D,L,xv,0,r,Ur,dUdr,F)
    for time in range(timesteps):
        
        if time > equitime:
            virialauto=0
            for i in range(N):
                virial += (1/2)*np.dot(dUdr[time,i,:],r[time,i,:])
                virialauto +=(1/2)*np.dot(dUdr[time,i,:],r[time,i,:])
            virialautocor[time-equitime-1]=virialauto

        xv[time + 1,:,:] = (xv[time,:,:] + vv[time,:,:]*h + ((h**2)/2)*F[time,:,:])%L
            
        for k in range(N):
            Ekin[time,:] = Ekin[time,:] + 0.5*m*(np.dot(vv[time,k,:],vv[time,k,:]))
        Energytotal[time,:] = Ur[time,:] + Ekin[time,:]
        
        F[time + 1,:,:] = force(N,D,L,xv,time+1,r,Ur,dUdr,F)
                                        
        vv[time + 1,:,:] = vv[time,:,:] + (h/2)*(F[time,:,:] + F[time+1,:,:])
               
        if (time % 10 == 0) & (time <= equitime):
            lam[time//10] = math.sqrt(Ekinequi/Ekin[time,:])
            vv[time+1,:,:] = lam[time//10] * vv[time+1,:,:]
    

    Tave= (2/3)*np.mean(Ekin[equitime+1:-1])/(N-1)
    Tautocor=(2/3)*Ekin[equitime+1:-1]/(N-1)
    Uave = np.sum(Ur[equitime:])/((timesteps-equitime+1)*N)
    
    kin2ave = np.mean(Ekin[equitime+1:-1]**2)
    kinave= np.mean(Ekin[equitime+1:-1])

    Cave = (3/2)*N/(1-(3/2)*N*((kin2ave-kinave**2)/(kinave**2)))

    compressibility = 1 - 1/(3*N*T) * (1/(timesteps-equitime+1)) * virial
    compressibilityautocor=1 - 1/(3*N*T) * virialautocor
    
    binwidth = round(L)/nbins
    
    for i in range(timesteps-equitime):  
        hist += np.histogram(r[equitime+1:,:,:][i][np.triu_indices(N,k=1)],bins = nbins, range = (0, round(L)))[0]
      

    hist = hist/(timesteps-equitime)
    rhist = np.linspace(binwidth/2,round(L)-binwidth/2,nbins)
    gr = ((2*L**3)/(N*(N-1)))*np.divide(hist,4*math.pi*binwidth*rhist**2)
 
    return Tave,Tautocor,Uave,Cave,compressibility,compressibilityautocor,r,Ur,dUdr,F,xv,vv,Ekin,Energytotal,rhist,hist,gr

