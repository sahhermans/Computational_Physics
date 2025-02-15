# Molecular dynamics simulation Argon - Computational Physics MSc Applied Physics 2020
# Authored by Hans Langerak & Sebastiaan Hermans 
# Date: 31/03/2020

# Main script
from IPython import get_ipython
get_ipython().magic('reset -sf')

from datetime import datetime
startTime = datetime.now()

import numpy as np
import matplotlib.pyplot as plt
import initialization as init 
import simulation as sim
import autocorrelation as aut
import plotfunction as pf
import datablocking as db

# Initialize all parameters of interest
D = 3                # number of dimensions
ncopies = 2          # number of copies of the unit cell of an FCC lattice in each dimension
N = 4*ncopies**D     # number of atoms 

timesteps = 5000     # number of timesteps of the simulation
equitime = 2500      # time untill equilibration is said to have been reach
h = 0.004            # length time step in unit time

T = 0.5                # dimensionless initial temperature
eoverkb = 119.8      # first fitting parameter Lennard-Jones potential for argon in [K]

sigma = 1.0          # second fitting parameter Lennard-Jones potential for argon in unit distance
rho = 1.2            # dimensionless density of argon in the simulation box 
L = (N/rho)**(1/3)   # size of the simulation box in unit distance
a = L/ncopies        # initialized lattice parameter of argon as defined by set parameters

m = 1.0              # mass argon atom in unit mass
nbins = round(L)*100 # number of bins for the pair correlation function

# Initialize positions & velocities
xv = init.fcc(N,D,a,ncopies,timesteps)
vv = init.maxwellboltzmann(T,N,D,timesteps)

Tave,Tautocor,Uave,Cave,compressibility,compressibilityautocor,r,Ur,dUdr,F,xv,vv,Ekin,Energytotal,rhist,hist,gr = sim.simulation(N,D,T,L,m,xv,vv,timesteps,equitime,h,nbins)

#Using datablocking to determine uncertainties
sdvcomp=db.datablocker(compressibilityautocor[49:-50])
sdvekin=db.datablocker(Ekin[equitime+1:-100])
sdvekin2=db.datablocker(Ekin[equitime+1:-100]**2)
sdvT=aut.autocorrelation(Tautocor)
sdvU=aut.autocorrelation(Ur[equitime+1:])

#Examples of figures that can be created using this script

#Plot of the energies 
#fig2=pf.plotfunction(range(timesteps),Energytotal[:-1],'Total energy',3,'t ($\sqrt{m\sigma^2 /\epsilon}$)','Energy ($\epsilon$)')
#fig3=pf.plotfunction(range(timesteps),Ur[:-1],'Potential energy', 3,'t ($\sqrt{m\sigma^2 /\epsilon}$)','Energy ($\epsilon$)')
#fig4=pf.plotfunction(range(timesteps),Ekin[:-1], 'Kinetic energy',3,'t ($\sqrt{m\sigma^2 /\epsilon}$)','Energy ($\epsilon$)')
#plt.savefig('392020energyT05report',dpi=200)

#Plot of the histogram
#fig5=pf.plotfunction(rhist,gr,'',5,'','')
#plt.savefig('392020grT05report',dpi=200)

#Boltzman distribution plot
#vplot = np.sort(np.sqrt(vv[0,:,0]**2 + vv[0,:,1]**2 + vv[0,:,2]**2))
#figboltzman=pf.histogram(vplot,15,'Initial velocity density plot','Absolute initial velocity [-]','Density [-]',6)
#plt.savefig('392020velocityhistT05report',dpi=200)

#ani=pf.animation3d(xv,L)
#ani

#Error Calculator C
ekinaverage=np.mean(Ekin[equitime:-1])
ekin2average=np.mean(Ekin[equitime:-1]**2)
delta=ekin2average-ekinaverage**2

covariance=np.mean((Ekin[equitime:]-np.mean(Ekin[equitime:])*(Ekin[equitime:]**2-np.mean(Ekin[equitime:]**2))))

udelta=np.sqrt(abs((sdvekin2)**2+(2*ekinaverage*sdvekin)**2-2*covariance))
udeltadk2=delta/ekin2average*np.sqrt((udelta/ekin2average)**2+(sdvekin**2/ekin2average**2)**2)

uc=Cave*udeltadk2*(2/(3*N)-udeltadk2)**(-2)/N
print(uc/N)
print(Cave/N-1.5)

print(datetime.now() - startTime)
