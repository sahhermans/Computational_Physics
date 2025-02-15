# Variational Quantum Monte Carlo - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 18/05/2020

# Main script

import numpy as np
import simulation as s
import functions as f

from datetime import datetime
startTime = datetime.now()

num_steps = 30000                   # Total number of walker steps.

walkers = 400                       # Total number of walkers.

displacement = 0.48                 # Value for walker displacement in metropolis 
                                    # algorithm.
                                    
equitime = 1                        # Number of steps after which we assume 
                                    # equilibrium position of the walkers.
                                    
alphas = np.array([0.15])           # Alpha values to be considered.

metropolis = s.metropolis_helium    # Choice for the metropolis algorithm to be 
                                    # used (s.metropolis_ + choice from helium, 
                                    # hydrogen & harmonic_oscillator).
                                    
dif = 0.001                         # Maximum difference between sebsequent alpha 
                                    # values in alpha minimization.
                                    
gamma = 0.5                         # Gamma value used for alpha minimization.



E,acceptanceratio,error_energies,var = s.energies(alphas,num_steps,walkers,displacement,
                                equitime,metropolis)

f.energies_plot_errorbar(alphas,E,error_energies)

print(datetime.now() - startTime)









