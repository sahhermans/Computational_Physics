# Variational Quantum Monte Carlo - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 21/05/2020

# Simulation script

import functions as f
import numpy as np
from numba import jit

@jit(nopython=True)
def metropolis_helium(num_steps,walkers,alpha,displacement,equitime):
    """
    Metropolis algorithm for the helium atom. 
    
     Input parameters:
        num_steps (int):            number of walker steps.
        
        walkers (int):              number of walkers.
        
        alpha (float):              alpha value considered.
        
        displacement (float):       value of the walker displacement.
        
        equitime (int):             number of steps after which we assume 
                                    equilibrium position of the walkers.
          
    Returned: 
        np.mean(E) (float):                 average energy of all walkers.
        
        np.mean(acceptanceratio) (float):   average acceptance ratio of walkers
                                            over total walk. 
        
        error_E (float):                    error in the energy
            
        var (float):                        variance in the energy
        
        np.mean(dyda) (float):              average derivative log of wavefunction 
                                            w.r.t. alpha.
        
        np.mean(eldyda) (float):            average of product of derivative and 
                                            the energy. 
    """
    
    E=np.zeros(walkers)
    E2=np.zeros(walkers)
    acceptanceratio=np.zeros(walkers)
    
    eldyda=np.zeros(walkers)
    dyda=np.zeros(walkers)
    
    r_trial1 = np.zeros(3)
    r_trial2 = np.zeros(3)
    
    for i in range(walkers):
        
        r1 = np.random.rand(3)*2 - 1
        r2 = np.random.rand(3)*2 - 1
        r12 = r1 - r2
        ov = np.array([0,0,0])
                    
        if np.all(r1 == ov) or np.all(r2 == ov) or np.all(r12 == ov): 
            raise ValueError("Initialized walker in origin")
    
        for j in range(num_steps):
            
            r_trial1 = r1 + displacement*(2*np.random.rand(3)-1)
            r_trial2 = r2 + displacement*(2*np.random.rand(3)-1)
            
            r_trial12 = r_trial1 - r_trial2
                        
            if (f.probability_den_helium(r_trial1,r_trial2,r_trial12,alpha)>=
               f.probability_den_helium(r1,r2,r12,alpha)):
                r1=r_trial1
                r2=r_trial2
                r12=r_trial12
                acceptanceratio[i] += 1/(num_steps)
                
            elif (f.probability_den_helium(r_trial1,r_trial2,r_trial12,alpha)<
                 f.probability_den_helium(r1,r2,r12,alpha)):
                if (f.probability_den_helium(r_trial1,r_trial2,r_trial12,alpha)/
                   f.probability_den_helium(r1,r2,r12,alpha) >= np.random.rand()):
                    r1=r_trial1
                    r2=r_trial2
                    r12=r_trial12
                    acceptanceratio[i] += 1/(num_steps)

            else:
                r1 = r1
                r2 = r2 
                r12 = r1-r2
            
            if j >= equitime:
                E[i] += (1/(num_steps-equitime))*(f.energy_loc_helium(r1,r2,r12,alpha))
                E2[i] += (1/(num_steps-equitime))*(f.energy_loc_helium(r1,r2,r12,alpha)**2)
      
                dyda[i] += (-(1/(num_steps-equitime))*(((np.linalg.norm(r12))**2)/2)
                           *(1+alpha*np.linalg.norm(r12))**(-2))
                eldyda[i] += (-(1/(num_steps-equitime))*f.energy_loc_helium(r1,r2,r12,alpha)*
                             (((np.linalg.norm(r12))**2)/2)*(1+alpha*np.linalg.norm(r12))
                             **(-2))

    var=np.mean(E2)-np.mean(E)**2
    error_E=np.sqrt(var)/walkers 
            
    return np.mean(E),np.mean(acceptanceratio),error_E,var,np.mean(dyda),np.mean(eldyda)

@jit(nopython=True)
def metropolis_hydrogen(num_steps,walkers,alpha,displacement,equitime):
    """
    Metropolis algorithm for the hydrogen atom.
    
    Input parameters:
        num_steps (int):            number of walker steps.
        
        walkers (int):              number of walkers.
        
        alpha (float):              alpha value considered.
        
        displacement (float):       value of the walker displacement.
        
        equitime (int):             number of steps after which we assume 
                                    equilibrium position of the walkers.
          
    Returned: 
        np.mean(E) (float):                 average energy of all walkers.
        
        np.mean(acceptanceratio) (float):   average acceptance ratio of walkers 
                                            over total walk. 
        
        error_E (float):                    error in the energy
            
        var (float):                        variance in the energy
        
        np.mean(dyda) (float):              average derivative log of wavefunction 
                                            w.r.t. alpha.
        
        np.mean(eldyda) (float):            average of product of derivative and 
                                            the energy. 
    """
    
    E=np.zeros(walkers)
    E2=np.zeros(walkers)
    acceptanceratio=np.zeros(walkers)
    eldyda=np.zeros(walkers)
    dyda=np.zeros(walkers)
    
    r_trial = np.zeros(3)
    
    for i in range(walkers):
        
        r = np.random.rand(3)*2 - 1
        ov = np.array([0,0,0])
        if np.all(r == ov): 
            raise ValueError("Initialized walker in origin")
        
        for j in range(num_steps):
            
            r_trial = r + displacement*(2*np.random.rand(3)-1)
            
            if f.probability_den_hydrogen(r_trial,alpha)>=f.probability_den_hydrogen(r,alpha):
                r=r_trial
                acceptanceratio[i] += 1/num_steps
    
            if f.probability_den_hydrogen(r_trial,alpha)<f.probability_den_hydrogen(r,alpha):
                if (f.probability_den_hydrogen(r_trial,alpha)/f.probability_den_hydrogen(r,alpha)
                   >=np.random.rand()):
                    r=r_trial
                    acceptanceratio[i] += 1/num_steps
        
            else:
                r=r
            
            if j >= equitime:
                E[i] += (1/(num_steps-equitime))*f.energy_loc_hydrogen(r,alpha)
                E2[i] += (1/(num_steps-equitime))*f.energy_loc_hydrogen(r,alpha)**2
                dyda[i] += -np.linalg.norm(r)/(num_steps-equitime)
                eldyda[i] += ((1/(num_steps-equitime))*f.energy_loc_hydrogen(r,alpha)
                             *(-np.linalg.norm(r)))
    
    var=np.mean(E2)-np.mean(E)**2
    error_E=np.sqrt(var)/walkers
       
    return np.mean(E),np.mean(acceptanceratio),error_E,var,np.mean(dyda),np.mean(eldyda)

@jit(nopython=True)
def metropolis_harmonic_oscillator(num_steps,walkers,alpha,displacement,equitime):
    """ 
    Metropolis algorithm for the harmonic oscillator.
    
    Input parameters:
        num_steps (int):            number of walker steps.
    
        walkers (int):              number of walkers.
        
        alpha (float):              alpha value considered.
        
        displacement (float):       value of the walker displacement.
        
        equitime (int):             number of steps after which we assume equilibrium 
                                    position of the walkers.
          
    Returned: 
        np.mean(E) (float):                 average energy of all walkers.
        
        np.mean(acceptanceratio) (float):   average acceptance ratio of walkers 
                                            over total walk. 
        
        error_E (float):                    error in the energy
            
        var (float):                        variance in the energy
        
        np.mean(dyda) (float):              average derivative log of wavefunction 
                                            w.r.t. alpha.
        
        np.mean(eldyda) (float):            average of product of derivative and 
                                            the energy. 
    """ 
    
    E=np.zeros(walkers)
    E2=np.zeros(walkers)
    acceptanceratio=np.zeros(walkers)
    eldyda=np.zeros(walkers)
    dyda=np.zeros(walkers)
      
    for i in range(walkers):
        
        x = np.random.rand()*2-1
        
        if x == 0: 
            raise ValueError("Initialized walker in origin")
        
        for j in range(num_steps):           
            
            x_trial = x + displacement*(2*np.random.rand()-1)
                        
            if (f.probability_den_harmonic_oscillator(x_trial,alpha)>
                f.probability_den_harmonic_oscillator(x,alpha)):
                x = x_trial
                acceptanceratio[i] += 1/num_steps
    
            elif (f.probability_den_harmonic_oscillator(x_trial,alpha)<
                  f.probability_den_harmonic_oscillator(x,alpha)):
                if (f.probability_den_harmonic_oscillator(x_trial,alpha)/
                    f.probability_den_harmonic_oscillator(x,alpha)
                   >= np.random.rand()):
                    x = x_trial
                    acceptanceratio[i] += 1/num_steps
        
            else:
                x=x
            
            if j >= equitime:
                E[i] += (1/(num_steps-equitime))*f.energy_loc_harmonic_oscillator(x,alpha)
                E2[i] += (1/(num_steps-equitime))*f.energy_loc_harmonic_oscillator(x,alpha)**2
                dyda[i] += -(x**2)/(num_steps-equitime)
                eldyda[i] += (1/(num_steps-equitime))*f.energy_loc_harmonic_oscillator(x,alpha)*(-x**2)
    
    var=np.mean(E2)-np.mean(E)**2
    error_E=np.sqrt(var)/walkers
            
    return np.mean(E),np.mean(acceptanceratio),error_E,var,np.mean(dyda),np.mean(eldyda)

def alpha_minimizer(num_steps,walkers,displacement,equitime,dif,gamma,alphas,metropolis):
    """ 
    Function used to find the minimum in energy as a function of parameter alpha. 
    
    Input parameters:
        num_steps (int):            number of walker steps.
        
        walkers (int):              number of walkers.
        
        alphas (float):             initial array of two alphas needed to start the minimization.
        
        displacement (float):       value of the walker displacement.
        
        equitime (int):             number of steps after which we assume equilibrium
                                    position of the walkers.
        
        dif (float):                maximum difference between sebsequent alpha 
                                    values in alpha minimization.
        
        gamma (float):              parameter used for alpha minimization.
        
        metropolis                  function of the metropolis algorithm you want to use,
                                    i.e. metropolis_harmonic_oscillator, metropolis_hydrogen
                                    or metropolis_helium
          
    Returned: 
        energies (float):           array of average energies at the considered 
                                    alphas.
        
        alphas (float):             array of alphas considered, last value gives 
                                    the energy minimizing one. 
        
        error_energies (float):     array of errors in average energies.
    """ 

    energies=np.array([]) 
    error_energies=np.array([])
    while np.abs(alphas[-1]-alphas[-2])>dif:
        E,acceptanceratio,error_E,var,dyda,eldyda=metropolis(num_steps,walkers,alphas[-1],
                                                      displacement,equitime)
        dEda=2*(eldyda-E*dyda)
        alpha_new=alphas[-1]-gamma*dEda
        alphas=np.append(alphas,np.array(alpha_new))
        energies=np.append(energies,np.array([E]))
        error_energies=np.append(error_energies,np.array([error_E]))
    
    alphas=np.delete(alphas,[alphas[0],alphas[1]])
    
    return alphas,energies,error_energies

    
@jit(nopython=True)
def energies(alphas,num_steps,walkers,displacement,equitime,metropolis):
    """ 
    Function used to ease calculations of the average energy for different values 
    of alpha.
    
    Input parameters:
        alphas (float):             alpha values considered.
        
        num_steps (int):            number of walker steps.
        
        walkers (int):              number of walkers.
        
        displacement (float):       value of the walker displacement.
        
        equitime (int):             number of steps after which we assume equilibrium 
                                    position of the walkers.
        
        dif (float):                maximum difference between sebsequent alpha 
                                    values in alpha minimization.
        
        gamma (float):              parameter used for alpha minimization.
          
    Returned: 
        energies (float):           array of average energies at the considered 
                                    alphas.
        
        alphas (float):             array of alphas considered, last value gives 
                                    the energy minimizing one. 
        
        error_E (float):            array of errors in average energies.
        
        var (float):                array of variances in average energies
    """ 

    energies=np.zeros(len(alphas))
    acceptanceratio=np.zeros(len(alphas))
    var=np.zeros(len(alphas))
    error_E=np.zeros(len(alphas))
    
    for i,alpha in enumerate(alphas):
        energies[i],acceptanceratio[i],error_E[i],var[i],dyda,eldyda=metropolis(num_steps,walkers,
                                                               alpha,displacement,
                                                               equitime)
    
    return energies,acceptanceratio,error_E,var

    
    
    
    
    