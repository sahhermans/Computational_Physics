# Computational Physics Project 3 - Hans Langerak & Sebastiaan Hermans

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import matplotlib.colors as col

from datetime import datetime
startTime = datetime.now()

L = 500                     # length 2D box
H = 100                     # height 2D box

deltat = 1                  # length time step
deltax = 1                  #length time step
u0 = np.ones((H,L))*0.15    # uniform initial velocity profile
#u0 = np.zeros((H,L))
#for p in range(L):         # Poiseuille initial velocity profile
#    u0[:,p] = 0.00001*np.arange(H)*(H-np.arange(H))
Re = 150
nu = 0.01 * L / Re
tau = 1/2 + (3*deltat*nu)/(deltax**2)
c = 1

e = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]]) # directional vectors of cube
w = (1/36) * np.array([16,4,1,4,1,4,1,4,1]) # weight directions
n = np.zeros((len(e),H,L))

for j in range(L):              # initialization profile
    for i in range(len(e)):
      n[i,:,j] = w[i]*(np.ones((H))+(3/(c**2))*np.dot(e[i][0],u0[:,0])+(9/(2*c**4))*(np.dot(e[i][0],u0[:,0]))**2-(3/(2*c**2)*(u0[:,0]*u0[:,0])))

boundary = np.zeros((H,L),dtype = bool) 
boundary[0,:] = True
boundary[-1,:] = True
boundary[int((H/3)), int(0):int(180)] = True
boundary[int((H/3)), int(260):int(L - 1)] = True
boundary[int((2*H/3)), int(0):int(180)] = True
boundary[int((2*H/3)), int(260):int(L - 1)] = True
boundary[int((0)):int((H/3 + 1)), int(180)] = True
boundary[int((2*H/3)):int(H), int(180)] = True

"""
#Circular boundary 
x = np.arange(0, 100)
y = np.arange(0, 500)
cx = 50  # Center of x-coordinate
cy = 30  # Center of y-coordinate
r = 10  # Radius of circle
circle = ((x[np.newaxis,:]-cx)**2 + (y[:,np.newaxis]-cy)**2 < r**2).transpose()
boundary[circle] = True
"""

bb = np.array([0,5,6,7,8,1,2,3,4])

# plot initial velocity profile
rho = sum(n)
ux = (n[1] + n[2] + n[8] - n[4] - n[5] - n[6])/rho
uy = (n[2] + n[3] + n[4] - n[6] - n[7] - n[8])/rho
ux[boundary] = 0
uy[boundary] = 0
plt.figure(figsize=(15, 5), dpi= 80)
plt.imshow(np.sqrt(ux**2+uy**2),origin='lower',cmap='coolwarm',interpolation='none')
plt.colorbar()
bImageArray = np.zeros((H, L, 4), np.uint8)
bImageArray[boundary,3] = 255								
plt.imshow(bImageArray, origin='lower',interpolation='none')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Initial velocity',fontsize = 14)

nn = np.zeros((len(e),H,L))
neq = np.zeros((len(e),H,L))
u = np.zeros((2,H,L))

theFig = plt.figure(figsize=(16,6))
ax = theFig.add_subplot(111,  xlim=(0, 500), ylim=(0, 100))
fluidImage = plt.imshow(np.sqrt((u[0])**2 + (u[1])**2), origin='lower',cmap='coolwarm', norm=plt.Normalize(0,0.25), interpolation='none')
bImageArray = np.zeros((H, L, 4), np.uint8)
bImageArray[boundary,3] = 255								
barrierImage = plt.imshow(bImageArray, origin='lower', interpolation='none')

#Time meter
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 1.1, '', transform=ax.transAxes)


def nextFrame(arg):			
    global boundary	
    for step in range(10):	   

        rho = sum(n)
        u[0] = (n[1] + n[2] + n[8] - n[4] - n[5] - n[6])/rho
        u[1] = (n[2] + n[3] + n[4] - n[6] - n[7] - n[8])/rho
        
        #if np.random.rand() > 0.99:                 # moving boundary
        #    boundary = np.roll(boundary,1,axis=1)
        
        for i in range(len(e)):
            # equilibrium solution densities
            neq[i] = (rho/tau)*w[i]*(1+(3/(c**2))*(e[i][0]*u[0] + e[i][1]*u[1])+(9/(2*c**4))*((e[i][0]*u[0]) + (e[i][1]*u[1]))**2-(3/(2*c**2)*((u[0])**2 + (u[1])**2)))
            
            # new density profile
            nn[i,:,:] = (1-1/tau)*n[i,:,:] + neq[i]
            nn[i,boundary] = n[bb[i],boundary]     
            
            # move densities
            n[i,:,:] = np.roll(np.roll(nn[i,:,:],e[i][0],axis=1),e[i][1],axis=0)
            
            # boundary outflow condition
            n[i,:,-1] = w[i]*(np.ones((H))+(3/(c**2))*np.dot(e[i][0],u0[:,0])+(9/(2*c**4))*(np.dot(e[i][0],u0[:,0]))**2-(3/(2*c**2)*(u0[:,0]*u0[:,0])))
            
    fluidImage.set_array(np.sqrt((u[0])**2 + (u[1])**2))
    
    time_text.set_text(time_template % (10*arg))

    return (fluidImage, barrierImage),time_text

animate = matplotlib.animation.FuncAnimation(theFig, nextFrame, interval=0)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar(fluidImage)
plt.show()
plt.title('Velocity animation',fontsize = 14)

def loop(iter):
    global u, u0, n, neq, nn, boundary, bb, tau, w, c, H, e
    for step in range(iter):	   
        
        rho = sum(n)
        u[0] = (n[1] + n[2] + n[8] - n[4] - n[5] - n[6])/rho
        u[1] = (n[2] + n[3] + n[4] - n[6] - n[7] - n[8])/rho
        
        #if np.random.rand() > 0.99:                 # moving boundary
        #    boundary = np.roll(boundary,1,axis=1)
        
        for i in range(len(e)):
            # equilibrium solution densities
            neq[i] = (rho/tau)*w[i]*(1+(3/(c**2))*(e[i][0]*u[0] + e[i][1]*u[1])+(9/(2*c**4))*((e[i][0]*u[0]) + (e[i][1]*u[1]))**2-(3/(2*c**2)*((u[0])**2 + (u[1])**2)))
            
            # new density profile
            nn[i,:,:] = (1-1/tau)*n[i,:,:] + neq[i]
            nn[i,boundary] = n[bb[i],boundary]     
            
            # move densities
            n[i,:,:] = np.roll(np.roll(nn[i,:,:],e[i][0],axis=1),e[i][1],axis=0)
            
            # boundary outflow condition
            n[i,:,-1] = w[i]*(np.ones((H))+(3/(c**2))*np.dot(e[i][0],u0[:,0])+(9/(2*c**4))*(np.dot(e[i][0],u0[:,0]))**2-(3/(2*c**2)*(u0[:,0]*u0[:,0])))
            
    return u

u = loop(7500)
plt.figure(figsize=(15, 5), dpi= 80)
Image = plt.imshow(np.sqrt(u[0]**2 +u[1]**2),cmap='coolwarm',origin='lower',interpolation='none',norm=plt.Normalize(0,np.max(np.sqrt(u[0]**2 +u[1]**2))))
bImageArray = np.zeros((H, L, 4), np.uint8)
bImageArray[boundary,3] = 255								
plt.imshow(bImageArray, origin='lower',interpolation='none')
plt.colorbar(Image)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Velocity plot',fontsize = 14)

u1 = np.mean(u[0,0:int(H/2),int(L/3)])
u2 = np.mean(u[0,0:int(H),int(L/2)])
A1 = int(H/2)
A2 = int(H) 
flow1 = u1*A1
flow2 = u2*A2

fig = plt.figure(figsize=(15,5), dpi= 80)
ax = plt.axes()
y,x = np.mgrid[0:H,0:L]
ustream = np.sqrt(u[0]**2 + u[1]**2)
lw = 2* ustream / np.max(ustream)
strm = ax.streamplot(x,y,u[0],u[1],color=np.sqrt(u[0]**2+u[1]**2),linewidth=lw,density=[5,1],cmap='coolwarm',norm = col.Normalize(0,np.max(np.sqrt(u[0]**2 +u[1]**2))))
plt.xlabel('x')
plt.ylabel('y')
bImageArray = np.zeros((H, L, 4), np.uint8)
bImageArray[boundary,3] = 255								
plt.imshow(bImageArray, origin='lower',interpolation='none')
fig.colorbar(strm.lines)
plt.title('Stream plot',fontsize = 14)

curl = np.roll(u[1],-1,axis=1) - np.roll(u[1],1,axis=1) - np.roll(u[0],-1,axis=0) + np.roll(u[0],1,axis=0)
plt.figure(figsize=(15, 5), dpi= 80)
Image = plt.imshow(curl,cmap='coolwarm',origin='lower',interpolation='none',norm=plt.Normalize(np.min(curl),np.max(curl)))
bImageArray = np.zeros((H, L, 4), np.uint8)
bImageArray[boundary,3] = 255								
plt.imshow(bImageArray, origin='lower',interpolation='none')
plt.colorbar(Image)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Curl plot',fontsize = 14)

print(datetime.now() - startTime)

