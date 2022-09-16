import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

##################################################################################
# OPIMIZATION USING TRF == TRUST REGION REFLECTIVE ALGO > ALLOWS BOUNDS          #
# HU & ZHOU 2014 > USE LM (Levenberg-Marquardt) WITHOUT BOUNDS FOR PRONY OPTIM   #
##################################################################################

# D(t) == CREEP COMPLIANCE

# D(t) = D0 + sum_i { D_i * ( 1 - exp(-t/10**rho_i) ) } + t*(1/nu_v)

def make_pronyJ(nn):
    def pronyJ(x,*p):
        prony=p[0]
        for i in range(1,2*nn,2):
            prony+=p[i]*(1-np.exp(-x/10**p[i+1]))
        prony+=x*p[2*nn+1]
        return prony
    return pronyJ

# E(t) == RELAXATION MODULUS (not DYNAMIC IT SEEMS...)         <-- JUSTIFIED BY INTERCONVERSION FROM CREEP COMPLIANCE MC

# E(t) = Eglassy + sum_i { E_i * ( exp(-t/10**rho_i) - 1 ) }
#      = Eglassy - sum_i { E_i * ( 1 - exp(-t/10**rho_i) ) }   <-- FORMULATION USED IN PSEUDO STRAIN/STRESS CALC

# Eglassy = Eequilibrium + sum_i { E_i }
#         = Erubbery     + sum_i { E_i }

# E(t) = Eequilibrium + sum_i { E_i * exp(-t/10**rho_i) }      <-- FORMULATION NEEDED FOR FINITE DIFF STRESS CALC

def make_pronyE(nn):
    def pronyE(x,*p):
        prony=p[0]
        for i in range(1,2*nn,2):
            prony+=p[i]*(np.exp(-x/10**p[i+1])-1)
        return prony
    return pronyE

#################################

def sigmoid(x,p):
    # > ttr + np.log10(JJ) + p0=[0.0]*4
    log_E=p[0]+p[1]/(1+np.exp(p[2]+p[3]*np.log10(1/x)))
    return log_E

def sigmoidLSQ(p,Xin,Eref):
    Esim = sigmoid(Xin,p)
    return np.sum((Esim-Eref)**2)/len(Eref)

def powerLaw(x,*p):
    # > p0=[0.0]*3
    power_D=p[0]+p[1]*x**p[2]
    return power_D

def MpowerLaw(x,*p):
    # > p0=[1.00]*4 + bounds=(0.0,np.inf) + dogbox
    power_D=p[0]+(p[2]-p[0])/((1+p[3]/x)**p[1])
    return power_D

def make_GpowerLaw(nn):
    # > p0=[1.0]*(2*nn+2) + bounds=(0.0,np.inf) + trf
    def GpowerLaw(x,*p):
        power_D=p[0]
        k=p[1]
        for i in range(2,2*nn+1,2):
            power_D+=p[i]/((1+p[i+1]/x)**k)
        return power_D
    return GpowerLaw
    
#################################

def reduced_t_order(tr,material):
    mat={'tr':tr,'J':material.iloc[:,2]}
    MC=pd.DataFrame(data=mat)
    MC.sort_values(by=['tr'], inplace=True)
    ttr=np.array(MC.iloc[:,0])
    MCarr=np.array(MC.iloc[:,1])

    return ttr,MCarr

def optim_prony(ttr,MCarr,typeI,nn):

    if typeI==1:
        fit_prony=make_pronyJ(nn)
        optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+2),bounds=(0,np.inf),method='trf')
    else:
        fit_prony=make_pronyE(nn)
        #optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+1),bounds=(0,np.inf),method='trf')
        optim_coeff,optim_stderr=curve_fit(fit_prony,ttr,MCarr,p0=[0.0]*(2*nn+1),method='trf')

    return optim_coeff,optim_stderr

def prony_plot(ttr,MCarr,typeI,coeff,stderr,nn):

    if typeI==1:
        print('\nProny Serie Parameters for Creep Compliance Master Curve = \n',coeff)
        fit_prony=make_pronyJ(nn)
        plt.title('Creep Compliance Master Curve vs. Reduced Time')
        plt.ylabel('Creep Compliance [1/GPa]')
    else:
        print('\nProny Serie Parameters for Erelax Master Curve = \n',coeff)
        fit_prony=make_pronyE(nn)
        plt.title('Relaxation Modulus Master Curve vs. Reduced Time')
        plt.ylabel('Relaxation Modulus [GPa]')
    
    mat_calc=fit_prony(ttr,*coeff)

    residuals = MCarr-mat_calc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((MCarr-np.mean(MCarr))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R² = '+str(r_squared))

    tsim=np.logspace(-8,8,50)
    mat_sim=fit_prony(tsim,*coeff)

    plt.plot(ttr,MCarr,'bs')
    plt.plot(ttr,mat_calc,'ro')
    plt.plot(tsim,mat_sim)
    plt.xlabel('Reduced Time [s]')
    plt.yscale('log')
    plt.xscale('log')
    #plt.grid(True)
    plt.show()