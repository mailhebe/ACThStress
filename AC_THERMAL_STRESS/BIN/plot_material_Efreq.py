import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_input(material,typeI,temp):
    
    if typeI==1:
        for k in temp:
            tempselection=material[material.iloc[:,0]==k]
            plt.plot(tempselection.iloc[:,1],tempselection.iloc[:,2])
        plt.title('Creep Compliance vs. Time')
        plt.xlabel('Time [s]')
        plt.ylabel('Creep Compliance [1/GPa]')
        plt.legend(temp,loc='best',title='Temperature [ºC]')
        #plt.grid(True)
        plt.show()
        
        for k in temp:
            tempselection=material[material.iloc[:,0]==k]
            plt.plot(np.log10(tempselection.iloc[:,1]),np.log10(tempselection.iloc[:,2]))
        plt.title('Log-Log Creep Compliance vs. Time')
        plt.xlabel('Log10(Time) [s]')
        plt.ylabel('Log10(Creep Compliance) [1/GPa]')
        plt.legend(temp,loc='best',title='Temperature [ºC]')
        #plt.grid(True)
        plt.show()

    else:
        
        for k in temp:
            tempselection=material[material.iloc[:,0]==k]
            plt.plot(tempselection.iloc[:,1],tempselection.iloc[:,2])
        plt.title('Dynamic Modulus vs. Frequency')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Dynamic Modulus [GPa]')
        plt.legend(temp,loc='best',title='Temperature [ºC]')
        #plt.grid(True)
        plt.show()

        for k in temp:
            tempselection=material[material.iloc[:,0]==k]
            plt.plot(np.log10(tempselection.iloc[:,1]),np.log10(tempselection.iloc[:,2]))
        plt.title('Log-Log Dynamic Modulus vs. Frequency')
        plt.xlabel('Log10(Frequency) [Hz]')
        plt.ylabel('Log10(Dynamic Modulus) [GPa]')
        plt.legend(temp,loc='best',title='Temperature [ºC]')
        #plt.grid(True)
        plt.show()