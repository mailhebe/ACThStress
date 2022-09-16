import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Williams-Landel-Ferry equation for reduced time calculation

def make_arrhenius(tref):
    def arrhenius(x,ea):
        return (ea/(2.303*8.314))*(1/(x+273.15)-1/(tref+273.15))
    return arrhenius

def optim_arrhenius(xmean,temp,tref,trefI):
    ylogat=np.insert(xmean,trefI,0)
    
    fit_arrhenius=make_arrhenius(tref)

    # Ea >=0 == Activatiopn Energy !
    #optim_coeff,optim_stderr=curve_fit(fit_arrhenius,temp,ylogat,bounds=(0,np.inf),method='trf')
    optim_coeff,optim_stderr=curve_fit(fit_arrhenius,temp,ylogat,method='trf')

    return optim_coeff,optim_stderr

def arrhenius_plot(temp,xmean,tref,trefI,coeff,stderr):
    ylogat=np.insert(xmean,trefI,0)

    fit_arrhenius=make_arrhenius(tref)
    logat_calc=fit_arrhenius(temp,*coeff)

    print('\nArrhenius model fit "log(aT)=Ea/(2.303R)*(1/T-1/Tref)" results :\n')
    print('Ea = ',coeff[0],'+/-',stderr[0,0]**0.5,'[J/mol]')

    residuals = ylogat- logat_calc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ylogat-np.mean(ylogat))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R² = '+str(r_squared))

    temprange=np.arange(-60,60,1)
    logat_sim=fit_arrhenius(temprange,*coeff)

    plt.plot(temp,ylogat,'bs')
    plt.plot(temprange,logat_sim,'g-')
    plt.title('Arrhenius Model - Parameter Optimization')
    plt.ylabel('Log10(aT) [s]')
    plt.xlabel('Temperature [ºC]')
    #plt.grid(True)
    plt.show()