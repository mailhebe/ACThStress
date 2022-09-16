import numpy as np
import matplotlib.pyplot as plt

def tref_sel(temp):

    tref=int(input('\nChose between one the following temperature as reference temperature Tref\n'+str(temp)+'\nSelection = '))
    trefInd=np.where(temp==tref)
    trefInd=int(trefInd[0])
    tempdel=np.delete(temp,trefInd)

    return tempdel,trefInd,tref

def mean_shift(material,typeI,tempdel,trefI,coeff):
    i=0
    xmean=np.zeros(len(tempdel))

    for k in tempdel:

        tempselection=material[material.iloc[:,0]==k]
        mat=np.log10(tempselection.iloc[:,2])
        
        if typeI==1:
            t=np.log10(tempselection.iloc[:,1])
        else:
            # test log(freq) optim for WLF parameters recovery #
            t=np.log10(tempselection.iloc[:,1])
            #t=np.log10(1/tempselection.iloc[:,1])
            ####################################################
            
        shift=(-coeff[trefI,1]+np.sqrt(coeff[trefI,1]**2-4*coeff[trefI,0]*(coeff[trefI,2]-mat)))/(2*coeff[trefI,0])
        xshift=t-shift

        xmean[i]=np.mean(xshift)

        i+=1
        
    return xmean

def equiv_slope(material,typeI,tempdel,tref,xmean):
    logtr=np.zeros(len(material.iloc[:,0]))
    tr=np.zeros(len(material.iloc[:,0]))

    for k in range(len(tr)):
        if material.iloc[k,0]==tref:
            if typeI==1:
                logtr[k]=np.log10(material.iloc[k,1])
            else:
                # test log(freq) optim for WLF parameters recovery #
                logtr[k]=np.log10(material.iloc[k,1])
                #logtr[k]=np.log10(1/material.iloc[k,1])
                ####################################################
        for tt in range(len(tempdel)):
            if material.iloc[k,0]==tempdel[tt]:
                if typeI==1:
                    logtr[k]=np.log10(material.iloc[k,1])-xmean[tt]
                else:
                    # test log(freq) optim for WLF parameters recovery #
                    logtr[k]=np.log10(material.iloc[k,1])-xmean[tt]
                    #logtr[k]=np.log10(1/material.iloc[k,1])-xmean[tt]
                    ####################################################

    tr=10**(logtr)
    mat=material.iloc[:,2]
    log_mat=np.log10(mat)

    plt.plot(tr,mat,'bo')

    if typeI==1:
        plt.title('Creep Compliance vs. Reduced Time')
        plt.ylabel('Log10(Creep Compliance) [1/GPa]')
        plt.xlabel('Reduced Time [s]')
    else:
        plt.title('Dynamic Modulus vs. Reduced Time')
        plt.ylabel('Log10(Dynamic Modulus) [GPa]')
        plt.xlabel('Reduced Frequency [Hz]')

    
    plt.yscale('log')
    plt.xscale('log')
    #plt.grid(True)
    plt.show()

    return tr