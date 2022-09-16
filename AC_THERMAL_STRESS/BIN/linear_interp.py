import numpy as np
import matplotlib.pyplot as plt

def timeinterp(time):

    # WE HAVE A PROBLEM WITH THE LINEAR INTERPOLATION, WHICH GENERATE 
    # ERRORS ON STRESS CALCULATION (DUE TO NUM FLUCTUATIONS?)
    # IF WE SKIP THE INTERPOLATION, WE HAVE GOOD RESULTS WITHOUT LOSS IN PRECISION!
    #tinterp=np.linspace(0,time[-1],(len(time)-1)*10+1)
    
    tinterp=time
    
    #tinterp[0]=time[0]
    #tinterp[0]=1/3600

    # NO NEED TO CORRECT t0 > WE USE DELTA_TIME / DELTA_REDTIME !

    tinterp[0]=0
    return tinterp

def tempinterp(time,tinterp,temp,depth):
    Tinterp=np.zeros((len(tinterp),len(depth)))
    time0=time
    time0[0]=0
    for k in range(len(depth)): Tinterp[:,k]=np.interp(tinterp,time0,temp[:,k])
    return Tinterp

def plotinterp(time,temp,tinterp,Tinterp):
    plt.plot(time[0:4],temp[0:4,0],'o')
    plt.plot(tinterp[0:31],Tinterp[0:31,0],'-x')
    plt.show()