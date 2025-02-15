# Project 2: Variational Quantum Monte Carlo

Authored by Hans Langerak (4592344) & Sebastiaan Hermans (4582608)

Here we present a variational monte carlo simulation for the helium atom.

For the scripts to work properly it is required that the user installs the Numba
and Numba-scipy packages.

On some computers the first time running the main file gives the value error:
ValueError: No function '_pyx_fuse_0pdtr' found in __pyx_capi_ of 'scipy.special.cython_special'

This problem is solved after running the script a second time. 

The repository contains multiple files:

`main.py` - Main script of the simulation. Here input parameters can be set, and desired output can be
    specified. Requires all other python scripts to run.

`functions.py` - Contains all the function required within the simulation file.

`simulation.py` - Performs the actual simulation using the input from the main script and using the functions 
    from the functions script.

`Report.pdf` - The final report on the results and performance of the presented simulation, made
    for the MSc course Computational Physics at Delft University of Technology.
    

Computational Physics, MSc Applied Physics 2020

Delft University of Technology

Date: 18/05/2020