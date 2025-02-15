# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on Wednesday** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented. 
2. write about your progress **before** the Tuesday in the next week with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
(due before 18 February)
Short plan for the week:
Firstly we used the lecture time to figure out how to use Gitlab and clone the created repository to our computers.
Secondly we made an introductory script in which we introduce all the required physics (constants, basic formulas etc.).
We subsequently want to meet on friday to start working on the milestones for next week. 
Beforehand we will have read the relevant chapter from dr. Thijssen's book.
We decided to not work individually, but only work together on the assignment for next week. 
Regarding code functionality, the aim is to implement all of the milestones described in week1.

Progress week 1:
On saturday we wrote a script (Script_1) in which all required functionalities for this week were implemented. 
In this first version we calculated the positions, velocities and the forces in separate x, y and z arrays 
for all timesteps. On monday we started working on a second script in which we stored the x, y and z 
components of the positions, velocities and forces in individual 3D matrices.

Firstly all parameters of interest are initialized. All parameters and their values are provided in 
SI units (m, kg etc.). 
Subsequently position and velocity arrays are initialized with numpy.zeros function with shape = (timesteps+1, number of particles).
Then the initial positions are initialized. For this a random number is picked with 
a random.rand function, which is then multiplied by L. This yield numbers between 0 and L (box size) for both x & y. 
The initial velocity is chosen with a random.normal function using 0 as mean as the thermal velocity (math.sqrt(kb*T/m))
as standard deviation. Next, by subtracting the mean of the generated vx and vy arrays from the same vx and vy arrays,
we make sure that the netto velocity of the particles within the box is zero.
After initializing arrays for r, Ur, dUdr and the force in the x & y directions,
the force on each particle is calculated the following way:

for time in range(timesteps):

    r[time,:] = np.sqrt((x[time,0]-x[time,1])**2 + (y[time,0]-y[time,1])**2) # distance between particles
    
    Ur[time,:] = 4*eoverkb*kb*((sigma/r[time,:])**12 - (sigma/r[time,:])**6) # Lennard-Jones potential
    
    dUdr[time,:] = 4*eoverkb*kb*( (-12*sigma**12)/(r[time,:]**13) + (6*sigma**6)/(r[time,:]**7) ) 

    Fx[time,0] = -dUdr[time,:] * (x[time,0]-x[time,1]) / r[time,:] # force on first particle in x-direction
    Fy[time,0] = -dUdr[time,:] * (y[time,0]-y[time,1]) / r[time,:] # force on first particle in y-direction
   
    Fx[time,1] = -dUdr[time,:] * (x[time,1]-x[time,0]) / r[time,:] # force on second particle in x-direction
    Fy[time,1] = -dUdr[time,:] * (y[time,1]-y[time,0]) / r[time,:] # force on second particle in y-direction

In the same for loop, a function for the total energy in the system in integrated:
Energytotal[time] = Ur[time,:] + (1/2)*m*(vx[time,0]**2 + vy[time,0]**2 + vx[time,1]**2 + vy[time,1]**2)

By means of a nested for-loop the Euler method is implented for time evolution:
    for particle in range(N):
                         
        x[time + 1,particle] = (x[time,particle] + vx[time,particle]*h) % L
        y[time + 1,particle] = (y[time,particle] + vy[time,particle]*h) % L
        #z[time + 1,particle] = (z[time,particle] + vz[time,particle]*h) % L
    
        vx[time + 1,particle] = vx[time,particle] + (1/m)*Fx[time,particle]*h
        vy[time + 1,particle] = vy[time,particle] + (1/m)*Fy[time,particle]*h
        #vz[time + 1,particle] = vz[time,particle] + (1/m)*Fz[time,particle]*h

Within this for-loop the periodic boundary condition was implemented by adding the % L after
the calculations for the positions in the next time step. We are not convinced that 
this is sufficient as this can mean that the calculated distance between two particles becomes larger
than should be the case. This is one of the questions we would like to ask this week. 
    
In the evening of monday I added Script_1.5: Within this new script, three different scripts can be found. 
The first contains the Script_1 script (separated in x, y and z components) only without the nested for-loop. 
The second script is an altered version of Script_1: now we work with 3D matrices containing the x, y and z components. 
This makes the script a bit faster. 
The third script contains the second script, only here we have started working with dimensionless numbers. 
For some reason it gives an 'int' not callable error for the Ur function -> need to find out why.

## Week 2
(due before 25 February)
Planning second week:
We will work on Wednesday 19/2/20 on the following problems:
Hans/Sebastiaan will derive the dimensionless kinetic energy individually and then compare.
We have already started implementing the dimensionless numbers, and we will finalize implementing them in the morning of wednesday.  
Next we will discuss how to implement the minimal image method, and try to implement it into the script. If time permits it we will also start working on the last milestone
if not then we will do that the following day.
We also decided that instead of one large script, from this week on we will start with a function script in which we define the required functions.

Progress on week 2:
On wednesday after the lecture we derived the expression for the dimensionless kinetic energy, we changed our existing molecular dynamics to now use dimensionless units 
and we implemented the minimal image convection and we wrote the code to make an animation of our simulation. On monday we figured out how to save
the animation and the figures. 
Note: all work has been done together.
The dimensionless kinetic energy is given by: Ek'=Ek/epsilon., however, we did not have to use this since we were already using dimensionless
parameters. The minimal image convection was implemented in the following way. We now use the function force to calculate the force:

def force(positionvector,time,N):
        
    for i in range(N):
        for j in range(N):
            if i != j:
                
                dv[time,i,:] = ((positionvector[time,i,:] - positionvector[time,j,:]) + L/2) % L - L/2
                                    
                r[time,:] = np.sqrt(np.dot(dv[time,i,:],dv[time,i,:]))
                Ur[time,:] = 4*((r[time,:])**(-12) - (r[time,:])**(-6)) 
                dUdr[time,:] = 4*(-12*r[time,:]**(-13) + 6*r[time,:]**(-7)) 

                F[time,i,:] = -dUdr[time,:] * dv[time,i,:] / r[time,:]

    return F[time,:,:]

In this function the minimal image convection is implented in the line for the positional difference vector dv[time,i,:]

In the following videos and figures we present our results:
First video:
Here we present an animation of two particles in 3D space with known initial positions velocities. We see exactly what we expected to see 
when the particles get close they repel eachother. Also the boundary condition seems to work because they never cross the boundary due to 
the repelling force. 
![alt text](/Videos/Plot_3D_towardseachother.mp4 "Animation of two particles in 3D space. Initial positions and velocities are chosen such that the particles move towards each other in 1D. We see what we expected: when the particles get too close, they start repelling each other such that they eventually move away from each other. The amount of timesteps is chosen so that the simulation does not yet break.")
This video was made with a limited amount of timesteps (1000). This is because after to many timesteps the scripts broke down due to a huge increase 
in velocity. We think that this is due to the euler method and will be fixed later on in the course when we move from this method. 
In the second video below you can see the simulation of two particles for random initial position and random intial velocities taken from a normal distribution
with zero mean and thermal velocity as standard deviation. What we see in the video is that the boundary conditions work as the particles emergy on
the other side of the volume. At 0:27 and a 1:05 you can find instances of where the particles get so close that they repel each other.
Looking at these videos can be quite confusing as its not clear when the particles are close to eachother as we only look at them from one angle.
![alt text](/Videos/Plot_3D_randomstartingpositions.mp4 "Animation of two particles in 3D space with random initial positions and velocities.")

And just for fun we include this video of our simulation for 30 particles:
![alt text](/Videos/Plot_3D_randomstartingpositions_30_particles.mp4 "Animation of 30 particles in 3D space with random initial positions and velocities.")


Finally, we also made some plots.
The first one below is a plot of the interatomic distance as a function of time in this case that the particles had known initial positions
and known initial velocities.
![alt text](/Images/interatomicdistance.png "Figure of the interatomic distance")
What we see in this figure is that there is clear period behavior in the first part which is in agreement with the first animation. However, after around a 
1000 timesteps the periodicity diminishes due to the huge increase in velocity due to the limitations of the euler method. Or due to the not rescaling of
the momenta of the particles

Furthermore, we also include of plot of the interatomic distance in the case of random initial conditions in the figure below:
![alt text](/Images/randominteratomicdistance.png "Figure of the interatomic distance")

Next is a figure where we plot the kinetic energy as a function of time for the same conditions as the plot above. We see that the kinetic energy
has a huge increase at the time we get errors due to the reason given above. 
![alt text](/Images/kineticenergy.png "Figure of the kinetic energy")
And here a plot of the kinetic energy for random initial conditions:
![alt text](/Images/randomkineticenergy.png "Figure of the kinetic energy")

Below is also a plot of the potential energy in the known case:
![alt text](/Images/potentialenergy.png "Figure of the potential energy")
Here is a plot of the potential energy in the case of random initial conditions.
![alt text](/Images/randompotentialenergy.png "Figure of the total energy")
We see clear delta peak like behavior when the particles get close to eachother in the case of random initial conditions which is in 
agreement with what we would expect based on the lennard jones potential.
Finally, a plot of the total energy is shown below. What we see is a huge increase in energy at the same time as the position get messy again for the same
reasons as given above. This is due to the contribution of the kinetic enregy.
![alt text](/Images/totalenergy.png "Figure of the total energy")
And here a plot of the total energy for random initial conditions:
![alt text](/Images/randomtotalenergy.png "Figure of the total energy")


There are multiple things we do not quite grasp yet. For exmaple we do not know why the kinetic energy increase just before they repel eachother. 
Also we would like to know if we can make this script work without our for loops. 

## Week 3
(due before 3 March)
Wednesday 26th of februari: We first are going to think about how we can eliminate the nested for loop that calculated the force. We think we can do 
this by making use of an array. Furthermore, our code does not recognize all the other particles when calculating the force. Therefore,
we will update the code using the array method given above to now support the force of all the particles. Important to consider is also that we 
figure out a way to prevent that the particles tries to calculate the force from itself on itself without making use of an if-statement.  
Finally, if time permits it today we will implement the velocity Verlett algorithm and make plots of the energies to verify if our simulation now
conserves energy a lot better then the Euler method.

Update 26th:
The nested for loops have been removed and we now use an array which results in a 300% increase in speed. Furthermore, we have also implented the 
Verlett algorithm and this has solved the problem that energy is not conserved. However, when we try to plot a lot more particles the code still 
manages to crash and energy is no longer conserved. The same happens for a large amount of time steps. We thought that this could be because the 
average velocity is not zero. We tried to correct for this by calculating the total mean velocity of all particles and subtracting it from all 
particles but it does not seem to work at this time. We would like to ask questions about other possibilities that could be responsible for this 
behaviour. Note, all this work has been done together to allow for more discussion about the right approach. 

Here is the code of the new way we implement the force calculation:
def force(positionvector,time,N):
        
    for i in range(N):
        dv = ((np.tile(xv[time,i,:],(N,1)) - xv[time,:,:]) + L/2) % L - L/2
        r[time,i,:] = np.sqrt(np.sum(dv**2,axis = 1))
        r[time,i,i] = 1
        Ur[time,:] = Ur[time,:] + np.sum((4*((r[time,i,:])**(-12) - (r[time,i,:])**(-6))))
        dUdr = 4*(-12*r[time,i,:]**(-13) + 6*r[time,i,:]**(-7)) 
        dUdr[i] = 0
        
        F[time,i,:] = np.dot(-dUdr/r[time,i,:],dv)

    return F[time,:,:]

We have removed the nested for loop by introducing an extra dimension in our difference vector array, given by dv.

Furthermore, here we present the implementation of the Verlett algorithm. Note that the small for loop at the end was intended to make sure the average 
velocity of the particles remained zero, however, this did not seem to work appropriately. For a very large amount of timesteps or for many particles in a tight box,
the simulation still went on to show way to large energies (maybe just don't simulate these boundary cases!). 

for time in range(timesteps):
    
    xv[time + 1,:,:] = (xv[time,:,:] + vv[time,:,:]*h + ((h**2)/2)*F[time,:,:])%L
        
    for k in range(N):
        Ekin[time,:] = Ekin[time,:] + 0.5*m*(np.dot(vv[time,k,:],vv[time,k,:]))
    Energytotal[time,:] = Ur[time,:]/2 + Ekin[time,:]
    
    F[time + 1,:,:] = force(xv,time + 1,N)
                                    
    vv[time + 1,:,:] = vv[time,:,:] + (h/2)*(F[time,:,:] + F[time+1,:,:])
    
    for p in range(D):
        vv[time+1,:,p] = vv[time+1,:,p] - stats.mean(vv[time+1,:,p])

As an illustration of the improvement below is a plot that contain the kinetic, potential and total energy for 6 particles. The first plot is for a known
situation and the second one for random initial conditions. It is clear that there is a huge improvement because now energy seems to be alot better conserved
from the straigt total energy line incontrast to the euler method. 

![alt text](/Images/totalenergyknownposvel.png "Figure of the energies for the new Verlett algorithm")

![alt text](/Images/totalenergyknownEuler.png "Figure of the energies in the case of the Euler approximation.")

Furtermore, we also present the same plot but now for random initial conditons.
![alt text](/Images/totalenergyrandom.png "Figure of the energies for the new Verlett algorithm.")

![alt text](/Images/totalenergyrandomEuler.png "Figure of the energies in the case of the Euler approximation.")


Finally, we also present animations for 6 particles using the improved code with known and random initial positions and velocities respectively..
![alt text](/Videos/plot_3D_scatter_animation_6particles_knownvelpos_v2.mp4 "Animation of six particles in 3D space. Initial positions and velocities are chosen such that the particles move towards each other in 1D. We see what we expected: when the particles get too close, they start repelling each other such that they eventually move away from each other. The amount of timesteps is chosen so that the simulation does not yet break.")

![alt text](/Videos/plot_3D_scatter_animation_6particles_randomposvel.mp4 "Animation of six particles in 3D space. Initial positions and velocities are chosen such that the particles move towards each other in 1D. We see what we expected: when the particles get too close, they start repelling each other such that they eventually move away from each other. The amount of timesteps is chosen so that the simulation does not yet break.")


## Week 4
(due before 10 March)
Planning this week:
We will work together untill around 14:00 on implementing the fcc latice as well as trying to show that the initial velocities obey a 
Maxwell-Boltzman distribution. Furthermore, if time permits we will perform the rescaling of temperature and show how the desired temperature 
is attained after a certain amount of rescaling and equilibrating. We will also start with studying the observables. The goal is to do at least
two of them this week. 

March 4th:
The fcc latice has been implemented using the code below, where ncopies is the amount of times the primitive unit cell is copied in each dimension.

```
N = 4*ncopies**D # number of atoms
L = (N/rho)**(1/3) # size of the simulation box in unit distance
def fcc(N,D,a,ncopies):
    """ 
    Returns an ffc lattice configuration. 
    
    Input parameters:
        N (int): gives the amount of particles in the system
        
        D (int): gives the amount of dimensions of the system
        
        a (float): gives the lattice constant of interest
        
        ncopies (int): gives the amount of copies of the unit cell 
                       of the FCC lattice in each dimension.
        
    Returned: 
        xv (float): array containing all x-, y- and z-components
        of all particles in the system.
    """
    xv = np.zeros((timesteps + 1, N, D), dtype = float)
    punitcell = np.array([[0,0,0],[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0]]) # primitive unit cell of FCC lattice
    particle = 0
    for g in range(len(punitcell)):
        for s in range(ncopies):
            for d in range(ncopies):
                for f in range(ncopies):
                    xv[0,particle,:] = punitcell[g] + a * np.array([s,d,f])
                    particle = particle + 1
                    
    return xv
```

Next, we show in the figure below that the initial velocities of the particles obey a maxwell boltzman distribution.
![alt text](/Images/392020velocityhist.png "Histogram of the absolute value of the velocity.")
This has been done using the following code:
```
T = 1 # temperature in [K]
vT = math.sqrt(T) # initial velocity in unit velocity
 
vv = np.zeros((timesteps + 1, N, D), dtype = float)
vv[0,:,:] = np.random.normal(0,vT,(N,D))
```

A figure in which we present the effect of rescaling the temperature is giving in the figure below:
![alt text](/Images/392020energy.png "Figure of the energies with rescaling. Total of 5000 timesteps of 0.004. Equilibrating every 10 timesteps untill equilibrium time 2500.")
The corresponding code for the rescaling is:
```
Ekinequi = (3/2)*(N-1)*T # kinetic energy given by equipartition theory

if (time % 10 == 0) & (time <= equitime):
        lam = math.sqrt(Ekinequi/Ekin[time,:])
        vv[time+1,:,:] = lam * vv[time+1,:,:]
```


We also calculated different observables using the following lines of code:

```
Tave = (2/3)*np.sum(Ekin[equitime:])/((N-1)*(timesteps-equitime))
Uave = np.sum(Ur[equitime:])/((timesteps-equitime+1)*N)

kin2ave = np.sum((Ekin[equitime:])**2)/(timesteps-equitime)
kinave = np.sum(Ekin[equitime:])/(timesteps-equitime)

Cave = (3/2)*N/(1-(3/2)*N*((kin2ave-kinave**2)/(kinave**2)))
BPrho = 1 - 1/(3*N*T) * (1/(timesteps-equitime+1)) * virial
```

with virial defined as
```
if time > equitime:
        for i in range(N):
            virial += (1/2)*np.dot(dUdr[time,i,:],r[time,i,:])
            
            
# pair correlation function
rave = np.ndarray.flatten(np.sum(r[equitime:,:,:],axis = 0)/(timesteps-equitime+1))
if np.sum(rave == 1.0) == N:
    rave = np.delete(rave,np.where(rave == 1))
else:
    print('Error: number 1 in r-array')

plt.figure(1)    
bw = 100
hist = plt.hist(rave,bins = round(L)*bw, range = (0,round(L)), weights = 1/2*np.ones(rave.shape))
binwidth = 1/bw
rhist = np.linspace(binwidth/2,round(L)-binwidth/2,round(L)*bw)
gr = ((2*L**3)/(N*(N-1)))*hist[0]/(4*math.pi*binwidth*rhist**2)
```

The following results have been obtained for an initial temperature of 1 and a density of 0.88.
The equilibrium temperature Tave is 0.99 which is in agreement with the example provided in jos thijsen book computational physics. 
The equilibrium potential energy is around -5.7 which is also in agreement with the literature in jos thijsens computational physics.
The specific heat is around 2.6. 
The kinetic energy is around 1.5 per particle.
The beta*pressure/rho flucates around 3.6.
Note that all these values vary per simulation and the given values are values we most often found. 
Whether the variance in these values is correct or not, we want to find out next week. We are especially unsure
about the found value for beta*pressure/rho, as this should be equal to 3 according to the book by dr. Thijssen.

Two results of the pair correlation are shown below for temperatures 0.5 and 1.
In the case of a temperature of 0.5 we would expect argon to be a solid. We can recognise
a solid in our result by the even spacing between peaks, as can be seen in the figure below (temperature 0.5):
![alt text](/Images/392020grT05.png "pair correlation function for a solid")
In the case of a liquid (or gas?) we see that there is less correlation, as can be seen in the figure below (temperature 1):
![alt text](/Images/392020grT1.png "pair correlation function for a liquid")

Finally, even tough we have implemented the pair correlation function we are not sure whether we 
did so correctly, as we binned the average of the distances between te particles over time. It could
also be that we should actually make a histogram at each timestep, and average these. We have seen that
this alters the result. We could not however manage to write this in such a way that python did not plot 
every histogram calculation, which made the plotting incredibly slow (python wanted to plot 2501 histograms).
Wednesday were going to ask questions about this. 

## Week 5
(due before 17 March)
This week we do distribute the tasks. Sebastiaan will focus on improving the calculation of the observables. While Hans will do the error calculation
for the observables. When these task are complete we will make a plan on how we want to validate our simulation. 

Edit: Finished the autocorellation function. I've added it as seperate notebook.Note that it needs to be cleaned up a bit but the script works.

The autocorrelation function has been implemented using the following code: 
```
def autocorellation(x):
    tijd=2498
    N=np.size(x)
    xa=np.zeros(tijd)
    for t in range(tijd):
        som=0
        for n in range(N-t):
            som=(x[n]*x[n+t])+som
        som1=np.sum(x[0:(N-t)])*np.sum(x[t:(N)])*(1/(N-t))
        xa[t]=(som-som1)*(1/(N-t))
    return xa/xa[0]
```

This code has been used to calculate the standard deviation for the heat capcacity and the pressure. Note, that we are going to make the autocorrelation
a bit more general since it now worked by implementing a number for tijd but this needs to be changed.

In the figure below we present the results of the autocorrelation fuctnion.
The heat capacity:
![alt text](/Images/correlationc.png "A fit that is used to extract the autocorrelation time which is used to derive the standard deviation")
Using the fit given above the autocorrelation time has been determined and used to find the standard deviation. The standarddeviation came out
to be: 0.1.

The pressure:
![alt text](/Images/correlationbprho.png "A fit that is used to extract the autocorrelation time which is used to derive the standard deviation")
Using the fit given above the autocorrelation time has been determined and used to find the standard deviation. The standarddeviation came out
to be: 0.25. Which is quite large given relative to the pressure. Maybe we determine the standard deviation in the wrong way for this problem,
we will ask questions about this later this week.

Make a plan for simulations to go into the report: How do you want to validate your simulation, and which observables/simulations do you want to run?

The plan for now is to simulate a box of 500 particles, and compare the results for the potential, kinetic and total energy per particle,
the final temperature of the system, the heat capacity and the pressure of the system and the pair correlation function. 
We want to do this for different phases of argon (different values of the initial temperature and density). 
Example initial temperatures and densities: Gas: T0 = 3, rho = 0.3, Liquid: T0 = 1, rho = 0.8, Solid: T0 = 0.5, rho = 1.2. 

We want to compare our results to the results found in the book by dr. Thijssen and in the Verlett paper, and possibly other papers found online. 
F.e. table 8.1 in the book by dr. Thijssen gives expected results for different values of rho for the potential energy per particle, the pressure 
of the system and the final temperature. We heard that the Verlett paper gives expected values for the heat capacity. The book by dr. Thijssen gives 
two different formulas for the heat capacity (page 210) with seemingly different dimensions. We don’t understand why, and would like to ask a TA on 
Wednesday what to do with this. Also, we are a bit confused by the pressure and heat capacity formulas, and don’t fully understand how to make them 
dimensionless with the factor kB. As our result for the pressure seems to overestimate the expected value this is definitely something we have to ask 
about this week.

We want to evaluate our result of the pair correlation function by looking up theory about the pair correlation function and comparing this with our results. 

At the moment we seem to be able to recreate the expected pair correlation functions for a gas, a liquid and a solid. We do however not understand 
why the peak of the most nearby particle has its mean slightly to the right of 1 sigma. As we are working with a Lennard-Jones potential we would 
expect that this mean is exactly at 1 sigma. Also, as we are only simulating one box we see a ‘cutoff’ distance after which the pair correlation 
function starts decreasing where you would expect it to continue. We have to find out this week whether it is expected of us to also simulate the 
26 boxes surrounding our box, so that we can fully recreate the theoretical pair correlation functions.

We have to think about also calculating the diffusion length, but as we already have a lot of observables that can be compared to theory we do not 
know if we are going to do so. 

For all observables mentioned above we have already implemented the calculation of the uncertainty, except for the pair correlation function, for 
which this calculation will be implemented somewhere this week. 

