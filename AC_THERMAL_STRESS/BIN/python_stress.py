import numpy as np
from math import exp

def stresscalc(tinterp,depth,Delta_RedTime,Etot,Delta_Etot,paramE,nn,Eref):

    Delta_maxwellstrain=np.zeros((nn,len(tinterp)-1,len(depth)))
    maxwellstrain=np.zeros((nn,len(tinterp),len(depth)))

    Delta_pseudostrain=np.zeros((len(tinterp)-1,len(depth)))
    pseudostrain=np.zeros((len(tinterp),len(depth)))

    stress=np.zeros((len(tinterp),len(depth)))

    for tt in range(0,len(tinterp)-1):

        for zz in range(len(depth)):
            Delta_pseudostrain[tt,zz]=paramE[0]*Delta_Etot[tt,zz]
            
            for k in range(1,2*nn,2):

                # "OverflowError: math range error" at tt=36
                # Check value of Delta_RedTime[tt,zz], since can cause overflow due to exponential function
                # YUP, inversion from positive to negative value of Delta_RedTime[36,0]≃-198294
                # If tweaked with 1/exp(+) ---> ZeroDivisionError: float division by zero

                #tempvar = Delta_RedTime[tt,zz]/(10**paramE[k+1])
                
                # ABSOLUTE VALUE GIVES QUITE GOOD VALUES ACTUALLY !
                tempvar = np.abs(Delta_RedTime[tt,zz]/(10**paramE[k+1]))
                
                Delta_maxwellstrain[int((k-1)/2),tt,zz]=(1-exp(-tempvar))*(Etot[tt,zz]-maxwellstrain[int((k-1)/2),tt,zz]) + \
                                                        1/tempvar*(tempvar+exp(-tempvar)-1)*Delta_Etot[tt,zz]

                maxwellstrain[int((k-1)/2),tt+1,zz]=maxwellstrain[int((k-1)/2),tt,zz]+Delta_maxwellstrain[int((k-1)/2),tt,zz]

                Delta_pseudostrain[tt,zz]-=paramE[k]*Delta_maxwellstrain[int((k-1)/2),tt,zz]
            
            pseudostrain[tt+1,zz]=pseudostrain[tt,zz]+Delta_pseudostrain[tt,zz]

            stress[tt+1,zz]=pseudostrain[tt+1,zz]*Eref

    return stress


def stresscalcFD(tinterp,depth,Delta_RedTime,Etot,Delta_Etot,paramE,nn):

    stress=np.zeros((len(tinterp),len(depth)))

    Eeq=paramE[0]
    for k in range(1,2*nn,2):
        Eeq-=paramE[k]

    for tt in range(0,len(tinterp)-1):

        for zz in range(len(depth)):
            stress[tt+1,zz]=Eeq*Etot[tt+1,zz]
            
            for k in range(1,2*nn,2):

                #stress[tt+1,zz]+=(paramE[k]*Delta_Etot[tt,zz]+stress[tt,zz])/(1+(Delta_RedTime[tt,zz]/10**paramE[k+1]))
                
                # ONCE AGAIN, CALCULATING WITH ABS VAL OF DELTA_REDTIME REDUCES ERRORS !
                stress[tt+1,zz]+=(paramE[k]*Delta_Etot[tt,zz]+stress[tt,zz])/(1+(np.abs(Delta_RedTime[tt,zz])/10**paramE[k+1]))

    return stress