import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def stressplot(tinterp,depth,stress):

    while True:

        view=0

        while (view<1) or (view>4):
            view=int(input('\nDo you want to see:\n1 - Stress VS Time at all depth?\n2 - Stress VS Time at given depth?\n3 - Stress profile at given time?'+ \
                            '\n4 - Punctual data : stress(time,depth)?\nSelection [1-4] = '))

        if view==1:
            tempstress=stress
            for k in range(len(depth)):
                plt.plot(tinterp,tempstress[:,k])
            plt.title('Calculated thermal stress VS Time')
            plt.xlabel('Time [h]')
            plt.ylabel('Thermal stress [GPa]')
            plt.legend(depth,loc='best',title='Depth [cm]')
            plt.show()
        elif view==2:
            print('\nDepth you want to plot?')
            for k in range(len(depth)):
                print(str(k)+' == '+str(depth[k])+'[cm]')
            zstress=int(input('\nSelection [0-'+str(k)+'] = '))
            tempstress=stress[:,zstress]
            plt.plot(tinterp,tempstress)
            plt.title('Calculated thermal stress VS Time')
            plt.xlabel('Time [h]')
            plt.ylabel('Thermal stress [GPa]')
            plt.show()
        elif view==3:
            tstress=int(input('\nTime you want to plot [h]?\n tstress = '))
            tindex=np.where(tstress>=tinterp)
            tempstress=stress[tindex[-1][-1],:]
            plt.plot(tempstress,-depth)
            plt.title('Calculated thermal stress VS Depth at t = '+str(tstress)+' [h]')
            plt.xlabel('Thermal stress [GPa]')
            plt.ylabel('Depth [cm]')
            plt.show()
            print('\nThermal stress at time t = '+str(tstress)+' [h] for all depth :\n')
            print(tempstress)
        else:
            tstress=int(input('\nTime you want to recover the data [h]?\n tstress = '))
            tindex=np.where(tstress>=tinterp)
            print('\nDepth you want to plot?')
            for k in range(len(depth)):
                print(str(k)+' == '+str(depth[k])+'[cm]')
            zstress=int(input('\nSelection [0-'+str(k)+'] = '))
            print('\n Thermal stress at time t = '+str(tstress)+' [h] and depth z = '+str(depth[zstress])+' [cm]')
            print('\n    ',str(stress[tindex[-1][-1],zstress]),' [GPa]')


        exportplot=str(input('\nDo you want to export the data? (y/N)\n'))

        if exportplot=='y' or exportplot=='Y':
            if view==1:
                exportftemp=pd.DataFrame(tempstress,index=[tinterp],columns=[depth])
                filepath = 'RESULTS/stressVStime_all.xlsx'
                exportftemp.to_excel(filepath,index=True,header=True)
            if view==2:
                exportftemp=pd.DataFrame(tempstress,index=[tinterp],columns=[zstress])
                filepath = 'RESULTS/stressVStime_'+str(depth[zstress])+'cm.xlsx'
                exportftemp.to_excel(filepath,index=True,header=True)
            if view==3:
                exportftemp=pd.DataFrame(tempstress,index=[tstress],columns=[depth])
                filepath = 'RESULTS/stressVSdepth_'+str(tstress)+'h.xlsx'
                exportftemp.to_excel(filepath,index=True,header=True)
        
        finished=str(input('\nFinished? (y/N)\n'))
        if finished=='y' or finished=='Y':break