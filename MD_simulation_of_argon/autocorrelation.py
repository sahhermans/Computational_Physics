import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def autocorrelation(x):
    
    """
    Calculates the standarddeviation of a given variable using autocorrelation. First the autocorrelation function
    is calculated, then the data is fitted to and exponential to get the autocorrelation time. This time is used
    to calculate the standarddeviation.
    
    Input parameters:
        x (array) : an array of variables 
            
    Returned:
        sigma (float)         : The standarddeviation of the input variable
        fig (figure)          : A figure presenting a fit of the autocorrelation data to a function exp(-t/tau)                   
        savedfigure (png)     : A saved figure of the fit
    """
    
    #Calculate Autocorrelation function
    time = 100
    N=np.size(x)
    xa=np.zeros(time)
    for t in range(time):
        summation=0
        for n in range(N-t):
            summation += (x[n]*x[n+t])
        summation1=np.sum(x[0:(N-t)])*np.sum(x[t:(N)])*(1/(N-t))
        xa[t]=(summation-summation1)*(1/(N-t))
    
    #Make a fit of the autocorrelation function
    def function(t,tau):
        return np.exp(-t/tau)
    
    t=np.arange(time)
    popt,pcov = curve_fit(function,t,xa/xa[0])
    
    A=np.mean(x**2)
    B=(np.mean(x))**2
    sigma=np.sqrt(2*popt[0]/np.size(xa)*(A-B))
    
    #This part can be activated to see the fit of the autocorrelation function to and exponent
    
    #fig=plt.figure()
    #plt.plot(t,np.exp(-t/popt[0]),label = 'fit ')
    #plt.plot(t,xa/xa[0],label = 'autocorrelation function')
    #plt.xlabel('t')
    #plt.ylabel('Xa')
    #plt.legend()
    #savedfigure=plt.savefig('correlationc',dpi=200)
    
    return sigma