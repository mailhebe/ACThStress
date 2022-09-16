import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Quadratic function to optimize on log(J)-log(t) data
def polylogJ(x,u,v,w):
    return u*x**2+v*x+w

# Linear function to optimize on log(E*)-log(t) data
def linlogE(x,u,v):
    return u*x+v

def optim_loglog(material,typeI,temp):

    i=0

    if typeI==1:
        optim_coeff=np.zeros((len(temp),3))
        optim_stderr=np.zeros((len(temp),3,3))
    else:
        optim_coeff=np.zeros((len(temp),2))
        optim_stderr=np.zeros((len(temp),2,2))

    for k in temp:

        tempselection=material[material.iloc[:,0]==k]
        mat=np.log10(tempselection.iloc[:,2])

        if typeI==1:
            t=np.log10(tempselection.iloc[:,1])
            coeff,pcov=curve_fit(polylogJ,t,mat)
        else:
            # test log(freq) optim for WLF parameters recovery #
            #t=np.log10(tempselection.iloc[:,1])
            t=np.log10(1/tempselection.iloc[:,1])
            ####################################################  

            coeff,pcov=curve_fit(linlogE,t,mat)

        optim_coeff[i,:]=coeff
        optim_stderr[i,:,:]=pcov
        
        i+=1

    return optim_coeff,optim_stderr

def optim_plot(material,typeI,temp,coeff,stderr):

    i=0
    temp2=np.zeros(2*len(temp))

    for k in temp:

        tempselection=material[material.iloc[:,0]==k]
        mat=np.log10(tempselection.iloc[:,2])

        if typeI==1:
            print('\nPolynomial fit "a.x²+b.x+c" for temperature =',k,'ºC\nResults :')
            print('a = ',coeff[i,0],'+/-',stderr[i,0,0]**0.5)
            print('b = ',coeff[i,1],'+/-',stderr[i,1,1]**0.5)
            print('c = ',coeff[i,2],'+/-',stderr[i,2,2]**0.5)
            t=np.log10(tempselection.iloc[:,1])
            mat_calc=polylogJ(t,*coeff[i,:])
        else:
            print('\nLinear fit "a.x+b" for temperature =',k,'ºC\nResults :')
            print('a = ',coeff[i,0],'+/-',stderr[i,0,0]**0.5)
            print('b = ',coeff[i,1],'+/-',stderr[i,1,1]**0.5)

            # test log(freq) optim for WLF parameters recovery #
            #t=np.log10(tempselection.iloc[:,1])
            t=np.log10(1/tempselection.iloc[:,1])
            ####################################################

            mat_calc=linlogE(t,*coeff[i,:])
                
        residuals=mat-mat_calc
        ss_res=np.sum(residuals**2)
        ss_tot=np.sum((mat-np.mean(mat))**2)
        r_squared=1-(ss_res/ss_tot)
        print('R² = '+str(r_squared))
        
        plt.plot(t,mat,'bs')
        plt.plot(t,mat_calc)

        temp2[2*i],temp2[2*i+1]=k,k

        i+=1

    if typeI==1:
        plt.title('Log-Log Creep Compliance vs. Time & Quadratic Optimization')
        plt.ylabel('Log10(Creep Compliance) [1/GPa]')
    else:
        plt.title('Log-Log Dynamic Modulus vs. Time & Linear Optimization')
        plt.ylabel('Log10(Dynamic Modulus) [GPa]')

    plt.xlabel('Log10(Time) [s]')
    plt.legend(temp2,loc='best',title='Temperature [ºC]')
    #plt.grid(True)
    plt.show()