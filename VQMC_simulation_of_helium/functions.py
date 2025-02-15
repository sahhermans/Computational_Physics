# Variational Quantum Monte Carlo - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 21/05/2020

# Functions script

import numpy as np
import matplotlib.pyplot as plt
from numba import jit

@jit(nopython=True)
def trialwf_helium(r1,r2,r12,alpha):
    """
    Generates a trial wavefunction for the helium atom.
    
    Input parameters:
        r1 (float):     current position vector of atom 1.
        
        r2 (float):     current position vector of atom 2.
        
        r12 (float):    current difference vector between positions of atoms 1 & 2.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        trial wave function of helium (float)
    """   
    
    return (np.exp(-2*np.linalg.norm(r1))*np.exp(-2*np.linalg.norm(r2)) *
           np.exp(np.linalg.norm(r12)/(2*(1+alpha*np.linalg.norm(r12)))))

@jit(nopython=True)
def energy_loc_helium(r1,r2,r12,alpha):
    """
    Computes the unitless local energy for the helium atom.
    
    Input parameters:
        r1 (float):     current position vector of atom 1.
        
        r2 (float):     current position vector of atom 2.
        
        r12 (float):    current difference vector between positions of atoms 1 & 2.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        local energy of helium (float)
    """
    El = (-4 + np.dot((r1/np.linalg.norm(r1) - r2/np.linalg.norm(r2)),(r1-r2)) * 
          1/(np.linalg.norm(r12)*(1+alpha*np.linalg.norm(r12))**2)-1/(np.linalg.norm(r12)
          *(1+alpha*np.linalg.norm(r12))**3)-1/(4*(1+alpha*np.linalg.norm(r12))**4)
          +1/np.linalg.norm(r12))
    return El

@jit(nopython=True)
def probability_den_helium(r1,r2,r12,alpha):
    """
    Computes the probability density of the trial wavefunction of the helium atom.
    
    Input parameters:
        r1 (float):     current position vector of atom 1.
        
        r2 (float):     current position vector of atom 2.
        
        r12 (float):    current difference vector between positions of atoms 1 & 2.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        probability density of the helium atoms (float)   
    """
    return trialwf_helium(r1,r2,r12,alpha)**2

@jit(nopython=True)
def trialwf_hydrogen(r,alpha):
    """
    Generates a trial wavefunction for the hydrogen atom.
    
    Input parameters:
        r (float):      current position vector of the atom.
           
        alpha (float):  alpha value considered.
        
    Returned: 
        trial wave function of hydrogen (float)
    """
    return np.exp(-alpha*np.linalg.norm(r))

@jit(nopython=True)
def energy_loc_hydrogen(r,alpha):
    """
    Computes the unitless local energy of the hydrogen atom.
    
    Input parameters:
        r (float):      current position vector of the atom.
           
        alpha (float):  alpha value considered.
        
    Returned: 
        local energy of hydrogen (float)
    """
    return -1/np.linalg.norm(r)-0.5*alpha*(alpha-2/np.linalg.norm(r))

@jit(nopython=True)
def probability_den_hydrogen(r,alpha):
    """
    Computes the probability density of the trial wavefunction of the hydrogen 
    atom.
    
    Input parameters:
        r (float):      current position vector of the atom.
           
        alpha (float):  alpha value considered.
        
    Returned: 
        probability density of hydrogen (float)
    """
    return trialwf_hydrogen(r,alpha)**2/(1/alpha)

@jit(nopython=True)
def trialwf_harmonic_oscillator(x,alpha):
    """
    Generates a trial wavefunction for the harmonic oscillator.
    
    Input parameters:
        x (float):      current position.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        trial wave function of the harmonic oscillator (float)
    """
    return np.exp(-alpha*x**2)

@jit(nopython=True)
def energy_loc_harmonic_oscillator(x,alpha):
    """
    Computes the unitless local energy of the harmonic oscillator.
    
    Input parameters:
        x (float):      current position.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        local energy of the harmonic oscillator (float)
    """
    return alpha+x**2*(0.5-2*alpha**2)

@jit(nopython=True)
def probability_den_harmonic_oscillator(x,alpha):
    """
    Computes the probability density of the trial wavefunction of the harmonic 
    oscillator.
    
    Input parameters:
        x (float):      current position.
        
        alpha (float):  alpha value considered.
        
    Returned: 
        probability density of the harmonic oscillator (float)
    """
    return trialwf_harmonic_oscillator(x,alpha)**2/np.sqrt(np.pi/(2*alpha))
    
def energies_plot_errorbar(alphas,E,error_energies):
    """
    Plots the computed energies as a function of alpha, with errorbars. 
    
    Input parameters:
        alphas (float): alpha values considered.
        
        E (float):      calculated energy values.
        
        stdE (float):   standard deviation in calculated energies
        
    Returned: 
        Plot of computed energy as a function of alpha, with errorbars
    """
    plt.errorbar(alphas,E,yerr=error_energies,marker='o',linestyle='', elinewidth=1,
                 ecolor='b',capsize=5)
    plt.xlabel(r'$\alpha$ (-)', fontsize= 16)
    plt.ylabel(r'$\langle E \rangle$ (-)', fontsize =16)













