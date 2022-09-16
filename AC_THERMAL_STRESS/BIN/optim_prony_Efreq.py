import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

def make_pronyJ(nn):
    def pronyJ(x,*p):
        prony=p[0]
        for i in range(1,2*nn,2):
            prony+=p[i]*(1-np.exp(-x/10**p[i+1]))
        prony+=x*p[2*nn+1]
        return prony
    return pronyJ

def make_pronyE(nn):
    def pronyE(x,*p):
        prony=p[0]
        for i in range(1,2*nn,2):
            prony+=p[i]*(np.exp(-x/10**p[i+1])-1)
        return prony
    return pronyE

def reduced_t_order(tr,material):
    mat={'tr':tr,'J':material.iloc[:,2]}
    MC=pd.DataFrame(data=mat)
    MC.sort_values(by=['tr'], inplace=True)
    ttr=np.array(MC.iloc[:,0])
    MCarr=np.array(MC.iloc[:,1])

    return ttr,MCarr

def optim_prony(ttr,MCarr,typeI):

    if typeI==1:
        nn=int(input('\nProny Serie terms for CCMC modeling?\n n = '))
        fit_prony=make_pronyJ(nn)
        optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+2),bounds=(0,np.inf),method='trf')
    else:
        nn=int(input('\nProny Serie terms for Erelax-MC modeling?\n n = '))
        fit_prony=make_pronyE(nn)
        #optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+1),bounds=(0,np.inf),method='trf')
        optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+1),method='trf')

    return optim_coeff,optim_stderr,nn

def prony_plot(ttr,MCarr,typeI,coeff,stderr,nn):

    if typeI==1:
        print('\nProny Serie Parameters for Creep Compliance Master Curve = \n',coeff)
        fit_prony=make_pronyJ(nn)
        plt.title('Creep Compliance Master Curve vs. Reduced Time')
        plt.ylabel('Creep Compliance [1/GPa]')
        
    else:
        print('\nProny Serie Parameters for Erelax Master Curve = \n',coeff)
        fit_prony=make_pronyE(nn)
        plt.title('Relaxation Modulus Master Curve vs. Reduced Frequency')
        plt.ylabel('Relaxation Modulus [GPa]')
    
    mat_calc=fit_prony(ttr,*coeff)
    
    residuals = MCarr-mat_calc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((MCarr-np.mean(MCarr))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R² = '+str(r_squared))

    tsim=np.logspace(-8,8,50)
    mat_sim=fit_prony(tsim,*coeff)

    if typeI==1:
        plt.plot(ttr,MCarr,'bs')
        plt.plot(ttr,mat_calc,'ro')
        plt.plot(tsim,mat_sim)
    else:
        plt.plot(1/ttr,MCarr,'bs')
        plt.plot(1/ttr,mat_calc,'ro')
        plt.plot(1/tsim,mat_sim)

    plt.xlabel('Reduced Time [s]')
    plt.yscale('log')
    plt.xscale('log')
    #plt.grid(True)
    plt.show()