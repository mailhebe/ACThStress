import math
import cmath
import numpy as np

def time_discr(t,dt):

    dhours=t*24
    tdiscr=int(dhours/dt+1)
    time=np.linspace(0,dhours,tdiscr,dtype=float)
    time[0]=1e-1 # added by Qingwen because of div/0 at first instant (at each depth)

    return time

def depth_discr(z,dz):

    if int(z)==z:
        length=z
    else:
        length=int(z)+1
        
    zsteps=np.arange(0,length+1,dz,dtype=float) # rounded depth steps in [cm]
    zsteps[-1]=z # real AC thickness value in [cm] at last element
    zm=zsteps/100 # depth steps in [m]

    return zm

def temp_profile_calculation(tvec,zvec,kac,kgb,aac,agb,tmax,tmin,tmean):

    p=np.array([12.8377+1.66606j,12.8377-1.66606j,12.2261+5.01272j,12.2261-5.01272j,
                10.9343+8.40967j,10.9343-8.40967j,8.77643+11.9219j,8.77643-11.9219j,
                5.22545+15.7295j,5.22545-15.7295j])

    m=np.array([-868.461+15457.4j, -868.461-15457.4j, 1551.63-8439.83j, 1551.63+8439.83j,
                -858.652+2322.07j, -858.652-2322.07j, 186.327-253.322j, 186.327+253.322j,
                -10.349+4.11094j, -10.349-4.11094j])

    ftemp=np.zeros((len(tvec),len(zvec)))

    for i in range(len(tvec)):
        for k in range(len(zvec)):
            ttemp=0
            for j in range(len(p)):
                temp0=p[j]/tvec[i] # here is the 0 div by time[0]
                
                r1=cmath.sqrt(temp0/kac)  # relaciona con cemento asfaltico
                r2=cmath.sqrt(temp0/kgb)  # relaciona con agregado
                mixprop=(aac*r1)/(agb*r2) # relaciona ambos
                
                ldiv=(1+mixprop.real)/(1-mixprop.real)
                
                fs=((tmax+tmin)/2)/temp0+(6*math.pi*(tmax-tmin))/(144*temp0**2+math.pi**2)
                a=(fs-tmean/temp0)/(1-ldiv*cmath.exp(2*zvec[-1]*r1))
                b=cmath.exp(zvec[k]*r1)-ldiv*cmath.exp((2*zvec[-1]-zvec[k])*r1)
                
                U=a*b
                FFs=temp0*U
                out=m[j]*FFs
                
                ttemp+=out.real
                
            ftemp[i,k]=ttemp+tmean

    return ftemp