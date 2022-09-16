import pandas as pd

from input import paramin, indexin

def input_material():
    
    datain1=''
    while datain1!=1 and datain1!=2:
        datain1=indexin('\nWhat do you want to import?\n1 - Creep Compliance Data\n2 - E* Data\nSelection = ')
        
    datain2=''
    while datain2!=1 and datain2!=2:
        datain2=indexin('\nWhat importation method do you want to use?\n1 - Excel Table Format\n2 - Clipboard Copy\nSelection = ')

    if datain1==1 and datain2==1:
        col_name=['temp_C','time_s','compliance_gpa^(-1)']
        datain=pd.read_excel('AC_THERMAL_STRESS/DATA/creep_compliance.xlsx')
        datain.columns=[col_name[0],col_name[1],col_name[2]]

    elif datain1==2 and datain2==1:
        col_name=['temp_C','freq_hz','Edyn_gpa']
        datain=pd.read_excel('AC_THERMAL_STRESS/DATA/dynamic_modulus.xlsx')
        datain.columns=[col_name[0],col_name[1],col_name[2]]

    elif datain1==1 and datain2==2:
        col_name=['temp_C','time_s','compliance_gpa^(-1)']
        ask='Please copy all creep compliance data {Temperature - Time - Compliance} from your Excel file!\nONCE your data is in the clipboard, press any key\n'
        entries='\n 1 - Temperature (C)\n 2 - Time (sec)\n 3 - Compliance (1/GPa)\n'
        datain=clipboard_input(col_name,ask,entries)

    else:
        col_name=['temp_C','freq_hz','Edyn_gpa']
        ask='Please copy all E* data {Temperature - Frequency - E*} from your Excel file!\nONCE your data is in the clipboard, press any key\n'
        entries='\n 1 - Temperature (C)\n 2 - Frequency (Hz)\n 3 - Dynamic Modulus (GPa)\n'
        datain=clipboard_input(col_name,ask,entries)

    return datain

def clipboard_input(col,ask,ent):

    while True:

        cask=input(ask)
        datain=pd.read_clipboard(header=None)

        if isinstance(datain.iloc[0,0], str):
            datain=datain.drop(datain.index[0])
            datain=datain.reset_index(drop=True)

        ctemp=indexin('The first column correspond to (give the number):'+ent)
        ctime_freq=indexin('The second column correspond to (give the number):'+ent)
        cmaterial=indexin('The third column correspond to (give the number):'+ent)

        datain.columns=[col[ctemp-1],col[ctime_freq-1],col[cmaterial-1]]
        datain=datain[[col[0],col[1],col[2]]]

        ### HERE CODE TO UNSURE WE'RE IMPORT FLOAT WITH POINT AS DECIMAL SEPARATOR
        ### BUT SEEMS USELESS, GENERALLY THAT'S BY DEFAULT...
        '''
        a=isinstance(datain.iloc[0,0],int)
        b=isinstance(datain.iloc[0,1],float)
        c=isinstance(datain.iloc[0,2],float)

        if a!=True:
            #datain[col[0]]=[x.replace(',','.') for x in datain[col[0]]]
            datain[col[0]]=datain[col[0]].astype(int)

        if b!=True:
            datain[col[1]]=[x.replace(',','.') for x in datain[col[1]]]
            datain[col[1]]=datain[col[1]].astype(float)

        if c!=True:
            datain[col[2]]=[x.replace(',','.') for x in datain[col[2]]]
            datain[col[2]]=datain[col[2]].astype(float)
        '''

        print('Please verify your data below.')
        print(datain)
        cask=input('Is is correct? (y/N)\n')
        if cask=='y' or cask=='Y': break

    return datain
            