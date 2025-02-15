# Readme Project 1: Molecular dynamics simulation of argon atoms using a 6-12 LJ potential
Authored by Hans Langerak (4592344) & Sebastiaan Hermans (4582608)

Here we present a molecular dynamics simulation of a system of up to 500 argon atoms, 
interacting through a 6-12 Lennard-Jones potential. The equations of motion are solved 
numerically using the velocity-Verlet algorithm as a symplectic integrator, ensuring energy 
conservation within our closed system. The minimal image convention is used to emulate an 
infinite system. 

The repository contains multiple files:

Main.py:
Main script of the simulation. Here input parameters can be set, and desired output can be
specified. Requires all other python scripts to run.

Initialization.py:
Initializes particles' positions and velocities in the way specified in the Main script.

Simulation.py:
Performs the actual simulation using the input from the main script and initialization
from the initialization script.

Autocorrelation.py:
Function that is used to calculate uncertainties in not time-averaged quantities.

Datablocking.py:
Function that is used to calculate uncertainties in time-averaged quantities.

Plotfunction.py:
Function used to plot calculated quantities and/or animate the simulated system.

Report.pdf:
The final report on the results and performance of the presented simulation, made 
for the MSc course Computation Physics at Delft University of Technology.

Journal.md:
The journal that was used to keep track of the progression made in the project. 

Videos (folder):
Contains 4 animations of the simulated system.


Computational Physics, MSc Applied Physics 2020

Delft University of Technology

Date: 01/04/2020









