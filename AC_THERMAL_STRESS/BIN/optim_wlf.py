import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Williams-Landel-Ferry equation for reduced time calculation

def make_wlf(tref):
    def wlf(x,c1,c2):
        return -c1*(x-tref)/(c2+x-tref)
    return wlf

def optim_wlf(xmean,temp,tref,trefI):
    ylogat=np.insert(xmean,trefI,0)
    
    fit_wlf=make_wlf(tref)

    # Normally, U and V are both >=0 in WLF model
    #optim_coeff,optim_stderr=curve_fit(fit_wlf,temp,ylogat,bounds=(0,np.inf),method='trf')
    # BUT bounding C1 and C2 value gives erroneous results...
    optim_coeff,optim_stderr=curve_fit(fit_wlf,temp,ylogat,method='trf')

    return optim_coeff,optim_stderr

def wlf_plot(temp,xmean,tref,trefI,coeff,stderr):
    ylogat=np.insert(xmean,trefI,0)

    fit_wlf=make_wlf(tref)
    logat_calc=fit_wlf(temp,*coeff)

    print('\nWLF model fit "log(aT)=-C1*(T-Tref)/(C2+(T-Tref))" results :\n')
    print('C1 = ',coeff[0],'+/-',stderr[0,0]**0.5)
    print('C2 = ',coeff[1],'+/-',stderr[1,1]**0.5,'[K]')

    residuals = ylogat- logat_calc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ylogat-np.mean(ylogat))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R² = '+str(r_squared))

    temprange=np.arange(-60,60,1)
    logat_sim=fit_wlf(temprange,*coeff)

    plt.plot(temp,ylogat,'bs')
    plt.plot(temprange,logat_sim,'g-')
    plt.title('Williams-Landel-Ferry Model - Parameter Optimization')
    plt.ylabel('Log10(aT) [s]')
    plt.xlabel('Temperature [ºC]')
    #plt.grid(True)
    plt.show()