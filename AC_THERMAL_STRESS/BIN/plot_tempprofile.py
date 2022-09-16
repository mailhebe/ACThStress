import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from input import paramin, indexin

def tempVStime(tempevol,timearray,deptharray):
    
    tstep=timearray[1]-timearray[0]

    while True:
        
        fullplot=str(input('\nDo you wanto to plot the whole graph (Y/n)\n'))

        startplot=0
        stopplot=timearray[-1]
        kplot=2

        if fullplot=='n' or fullplot=='N':
            halfplot=str(input('\nDo you wanto to plot the graph once equilibrium is reached? (Y/n)\n'))
            kplot=1
            if halfplot=='n' or halfplot=='N':
                startplot=indexin('\nChoose a starting time [h] :')
                stopplot=indexin('\nChoose a ending time [h] :')
                kplot=0

        if kplot==2:
            time2=timearray
            time2[0]=0
            ftemp2=tempevol
        elif kplot==1:
            time2=timearray[int(24/tstep):]-24 # Recover only second day and correct time stamp
            ftemp2=tempevol[int(24/tstep):,:] # Avoid discontinuity due to div/0
        else:
            time2=timearray[int(startplot/tstep):int(stopplot/tstep)+1]
            ftemp2=tempevol[int(startplot/tstep):int(stopplot/tstep)+1,:]
        
        plt.plot(time2,ftemp2)
        plt.title('Temperature evolution for different depth')
        plt.xlabel('Time [h]')
        plt.xticks(np.arange(time2[0],time2[-1]+1,tstep*2))
        plt.xlim(time2[0],time2[-1])
        plt.ylabel('Temperature [ºC]')
        plt.ylim(round(min(ftemp2[:,0]))-1,round(max(ftemp2[:,0]))+1)
        #plt.legend(deptharray,loc='center right',bbox_to_anchor=(1.25,0.5),title='Depth [cm]')
        plt.legend(deptharray,loc='best',title='Depth [cm]')
        plt.grid(True)
        plt.show()

        exportplot=str(input('\nDo you want to export the data? (y/N)\n'))

        if exportplot=='y' or exportplot=='Y':
            exportftemp=pd.DataFrame(ftemp2,index=[time2],columns=[deptharray])
            filepath = 'RESULTS/temp_profile.xlsx'
            exportftemp.to_excel(filepath,index=True,header=True)
        
        finished=str(input('\nFinished? (y/N)\n'))
        if finished=='y' or finished=='Y':break

def depthVStemp(tempevol,tstep,deptharray):
    
    while True:

        ttime=indexin('\nHour of the temperature profile wanted : ')

        plt.plot(tempevol[int(ttime/tstep),:],-deptharray) # Depth vs Temp
        plt.ylabel('Depth [cm]')
        #plt.ylim(round(-deptharray[-1])-0.5,0.5)
        #plt.plot(tempevol[int(ttime/tstep),:],-deptharray*0.3937008) # Depth_inches vs Temp
        #plt.ylabel('Depth [in]')
        #plt.ylim(round(-deptharray[-1]*0.3937008)-0.5,0.5)

        plt.title('Temperature profile at '+str(ttime)+' h')
        plt.xlabel('Temperature [ºC]')
        plt.grid(True)
        plt.show()

        # export in excel with automatic naming // hour

        exportplot=str(input('\nDo you want to export the data? (y/N)\n'))

        if exportplot=='y' or exportplot=='Y':
            exportftemp=pd.DataFrame(tempevol[int(ttime/tstep),:],index=[deptharray],columns=[ttime])
            filepath = 'RESULTS/depthVStemp_'+str(ttime)+'h.xlsx'
            exportftemp.to_excel(filepath,index=True,header=True)
        
        finished=str(input('\nFinished? (y/N)\n'))
        if finished=='y' or finished=='Y':break