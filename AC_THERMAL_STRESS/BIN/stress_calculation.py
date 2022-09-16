import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from optim_arrh import make_arrhenius
from optim_wlf import make_wlf

def tsecWLF(Tinterp,coeff,tref):
    func_wlf=make_wlf(tref)
    lat=func_wlf(Tinterp,*coeff)
    return lat

def tsecArrh(Tinterp,coeff,tref):
    func_arrhenius=make_arrhenius(tref)
    lat=func_arrhenius(Tinterp,*coeff)
    return lat

def stressprecalc(tinterp,Tinterp,depth,lat,bmix,tmean):

    tsec=tinterp*3600

    RedTime=np.zeros((len(tinterp),len(depth)))
    Delta_RedTime=np.zeros((len(tinterp)-1,len(depth)))
    Etot=np.zeros((len(tinterp),len(depth)))
    Delta_Etot=np.zeros((len(tinterp)-1,len(depth)))

    for k in range(len(depth)):   
        RedTime[:,k]=tsec/10**(lat[:,k])            # == reduced_time
        Delta_RedTime[:,k]=np.diff(RedTime[:,k])    # == Delta_redtime
        Etot[:,k]=bmix*(Tinterp[:,k]-tmean)         # == total_strain = alpha * [ T(t,z) - Tequi ]
        Etot[0,k]=0                                 # CORRECTION ERROR OF TH. PROFILE :
                                                    #   AT t=0 > SLAB AT EQUILIBRIUM TEMP (T=Tmean)
                                                    #   THEN Etot(t=0) = 0
                                                    #   SOLVE OFFSET AT START OF PSEUDO STRAIN CALC
        Delta_Etot[:,k]=np.diff(Etot[:,k])          # == Delta_total_strain
    
    '''
    exportftemp=pd.DataFrame(lat,index=[tinterp],columns=[depth])
    filepath = 'RESULTS/lat.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)
    '''

    return RedTime,Delta_RedTime,Etot,Delta_Etot