import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from optim_prony import make_pronyJ,make_pronyE
from optim_prony import sigmoid,powerLaw,MpowerLaw
from optim_prony import make_GpowerLaw

def make_dpronyJ(nn):
    def dpronyJ(x,*p):
        dprony=0
        for i in range(1,2*nn+1,2):
            dprony+=(p[i]/10**p[i+1])*np.exp(-x/10**p[i+1])
        dprony+=p[2*nn+1]
        return dprony
    return dpronyJ

def interconv_CCMC_Edyn(ttr,MCarr,coeff,nn):

    dpronyJ1=make_dpronyJ(nn)
    dpronyJ=dpronyJ1(ttr,*coeff)
    
    fit_prony=make_pronyJ(nn)
    Jcalc=fit_prony(ttr,*coeff)

    n=np.abs(ttr*dpronyJ/(Jcalc))
    alpha=(np.sin(n*np.pi)/(n*np.pi))**(1/n)
    Dtalpha=fit_prony(ttr/alpha,*coeff)
    E=1/Dtalpha

    return E

def general_interconv_prony(ttr,MCarr,coeff,nn):

    fit_prony=make_pronyJ(nn)
    Jcalc=fit_prony(ttr,*coeff)

    ttr2=ttr[0:len(ttr)-1]
    nlog=np.zeros(len(Jcalc)-1)

    for k in range(len(Jcalc)-1):
        nlog[k]=np.abs((np.log10(Jcalc[k])-np.log10(Jcalc[k+1]))/(np.log10(ttr[k])-np.log10(ttr[k+1])))
    
    alphalog=(np.sin(nlog*np.pi)/(nlog*np.pi))**(1/nlog)
    Dtalphalog=fit_prony(ttr2/alphalog,*coeff)
    Elog=1/Dtalphalog

    return Elog

def general_interconv(tsim,MCsim,coeff,nn,func):

    ttemp=tsim[0:len(tsim)-1]
    nlog=np.zeros(len(MCsim)-1)

    for k in range(len(MCsim)-1):
        if func=='sigmoid':
            nlog[k]=np.abs((np.log10(10**MCsim[k])-np.log10(10**MCsim[k+1]))/(np.log10(tsim[k])-np.log10(tsim[k+1])))
        else:
            nlog[k]=np.abs((np.log10(MCsim[k])-np.log10(MCsim[k+1]))/(np.log10(tsim[k])-np.log10(tsim[k+1])))
    
    alphalog=(np.sin(nlog*np.pi)/(nlog*np.pi))**(1/nlog)

    if func=='prony':
        fit_prony=make_pronyJ(nn)
        Dtalphalog=fit_prony(ttemp/alphalog,*coeff)
    elif func=='sigmoid':
        Dtalphalog=sigmoid(ttemp/alphalog,*coeff)
    elif func=='power':
        Dtalphalog=powerLaw(ttemp/alphalog,*coeff)
    elif func=='Mpower':
        Dtalphalog=MpowerLaw(ttemp/alphalog,*coeff)
    elif func=='GMpower':
        fit_GMPL=make_GpowerLaw(nn)
        Dtalphalog=fit_GMPL(ttemp/alphalog,*coeff)

    if func=='sigmoid':
        Elog=1/(10**Dtalphalog)
    else:
        Elog=1/Dtalphalog    

    return Elog

def optim_interconv(ttr,E,nn):

    #nn=int(input('\nProny Serie terms for E*-MC modeling?\n n = '))

    fit_pronyE=make_pronyE(nn)
    optim_coeff,optim_stderr=curve_fit(fit_pronyE,ttr,E,p0=[0.0]*(2*nn+1),bounds=(0,np.inf),method='trf')
    
    return optim_coeff,optim_stderr


def interconv_plot(ttr,E,coeff,stderr,nn):

    print('\nProny Serie Parameters for Interconverted E* Master Curve = \n',coeff)

    fit_pronyE=make_pronyE(nn)
    Ecalc=fit_pronyE(ttr,*coeff)

    residuals = E-Ecalc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((E-np.mean(E))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R² = '+str(r_squared))

    tsim=np.logspace(-4,7,50)
    Esim=fit_pronyE(tsim,*coeff)

    plt.plot(ttr,E,'bs')
    plt.plot(ttr,Ecalc,'ro')
    plt.plot(tsim,Esim,'g')
    plt.title('Dynamic Modulus Master Curve vs. Reduced Time')
    plt.ylabel('Dynamic Modulus [GPa]')
    plt.xlabel('Reduced Time [s]')
    plt.yscale('log')
    plt.xscale('log')
    #plt.grid(True)
    plt.show()