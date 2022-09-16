# # # # # # # # # # # # # # # # # # # # # # # #
#                                             # 
#     ACTS Calc: Asphalt Concrete Thermal     #
#                   Stress Calculator         #
#            By: Benjamin J. Mailhé           #
#                 21-12-2021                  #
#            Version Beta-release             #  
#                                             #
# # # # # # # # # # # # # # # # # # # # # # # #

######################################
#   IMPORTATION OF USEFUL PACKAGES   #
######################################
import math
import cmath
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit,least_squares
# from scipy.optimize import minimize
# import matplotlib.pyplot as plt
# from matplotlib.pyplot import figure

##################################
#   IMPORTATION OF GUI LIBRARY   #
##################################

####### dearpygui v0.6.415 #######
from dearpygui.core import *
from dearpygui.simple import *

####### for latest version #######
# from dearpygui.dearpygui import *

###########################
#   LINK TO SUBROUTINES   #
###########################
import sys
import os
import webbrowser

sys.path.insert(1,'AC_THERMAL_STRESS/BIN/')
path=os.path.abspath('')

############################
#   FUNCTION IMPORTATION   #
############################
from input import paramin, indexin
from temp_profile_gui2007 import wind_approx, wind_poly_val, r2_calc
from temp_profile_gui2007 import conversion_sewing
from temp_profile_gui2007 import pavement_depth_const
from temp_profile_gui2007 import thermal_profile_calc
from optim_rawmat import optim_loglog,polylogJ,linlogE
from tref_shift import mean_shift,equiv_slope
from optim_prony import reduced_t_order
from optim_wlf import optim_wlf,make_wlf
from optim_arrh import optim_arrhenius,make_arrhenius
from optim_prony import make_pronyJ,make_pronyE
from optim_prony import optim_prony
from optim_prony import sigmoid,sigmoidLSQ,powerLaw,MpowerLaw
from optim_prony import make_GpowerLaw
from interconvertion import interconv_CCMC_Edyn,optim_interconv
from interconvertion import general_interconv
from linear_interp import timeinterp,tempinterp
from stress_calculation import tsecWLF,tsecArrh,stressprecalc
from python_stress import stresscalc, stresscalcFD

#############################################################################
# CONSTANTS
#############################################################################

stefan_boltzmann=5.67e-8 # [W/m²ºK⁴]

ref_temp_flow=[300,350]
ref_viscosity_flow=[1.589e-05,2.092e-05]
ref_conductivity_flow=[2.630e-02,3.000e-02]
ref_diffusivity_flow=[2.250e-05,2.990e-05]
ref_prandtl_flow=[0.707,0.700]

#############################################################################
def start(sender, data):
    show_item('Weather Data Input')
    hide_item('Start Page')

def about(sender, data):
    # show_item('About')
    # hide_item('Start Page')
    webbrowser.open_new(path+'/AC_THERMAL_STRESS/INTER_DATA/About.pdf')

def readme(sender, data):
    webbrowser.open_new(path+'/AC_THERMAL_STRESS/INTER_DATA/Readme_COBMok.pdf')

def start2(sender, data):
    show_item('Start Page')
    hide_item('Weather Data Input')
    hide_item('About')

#############################################################################
def input_weather_data(sender, data):
    days_to_compute=int(get_value('days_compute'))

    if days_to_compute==1:
        configure_item('day_data##1',items=['Day 1'],show=True)
        configure_item('day_data##2',items=['Day 1'])
    elif days_to_compute==2:
        configure_item('day_data##1',items=['Day 1','Day 2'],show=True)
        configure_item('day_data##2',items=['Day 1','Day 2'])
    elif days_to_compute==3:
        configure_item('day_data##1',items=['Day 1','Day 2','Day 3'],show=True)
        configure_item('day_data##2',items=['Day 1','Day 2','Day 3'])

def weather_test(sender, data):

    days_to_compute=int(get_value('days_compute'))

    days=[*range(0,days_to_compute)]
    weather_data=pd.read_excel(path+'/DATA/weather.xlsx',sheet_name=days,header=0)

    if data=='Day 1':
        table_in=weather_data[0].values
        table_in=np.round(table_in,2)
        table_in=table_in.tolist()
    elif data=='Day 2':
        table_in=weather_data[1].values
        table_in=np.round(table_in,2)
        table_in=table_in.tolist()
    elif data=='Day 3':
        table_in=weather_data[2].values
        table_in=np.round(table_in,2)
        table_in=table_in.tolist()

    clear_table('weather_table')
    set_headers('weather_table',weather_data[0].columns.tolist())
    set_table_data('weather_table',table_in)
    show_item('weather_table')

    show_item('Next##1')

def show(sender, data):
    
    days_to_compute=int(get_value('days_compute'))
    days=[*range(0,days_to_compute)]
    weather_data=pd.read_excel(path+'/DATA/weather.xlsx',sheet_name=days,header=0)

    wind_param=np.zeros((days_to_compute,4))
    wind_stderr=np.zeros((days_to_compute,4,4))
    wind_param,wind_stderr=wind_approx(days_to_compute,weather_data)

    day_col=[]
    for i in range(days_to_compute): day_col.append('Day '+str(i+1))

    wind_poly=wind_poly_val(days_to_compute,weather_data,wind_param)
    exportftemp=pd.DataFrame(wind_poly.transpose(),columns=day_col)
    filepath = path+'/RESULTS/wind_poly.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    wind_r2=r2_calc(days_to_compute,weather_data,wind_poly,wind_stderr)
    wind_export=np.c_[wind_param,wind_r2]
    exportftemp=pd.DataFrame(wind_export,columns=['a','b','c','d','R²'],index=day_col)
    filepath = path+'/RESULTS/wind_param.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    hide_item('Weather Data Input')
    show_item('Wind Speed Approximation')

#############################################################################
def wind_speed_plot(sender, data):

    wind_poly=pd.read_excel(path+'/RESULTS/wind_poly.xlsx',header=0,index_col=0)
    wind_param=pd.read_excel(path+'/RESULTS/wind_param.xlsx',header=0,index_col=0)

    days_to_compute=int(get_value('days_compute'))
    days=[*range(0,days_to_compute)]
    weather_data=pd.read_excel(path+'/DATA/weather.xlsx',sheet_name=days,header=0)
    
    if data=='Day 1':
        data_raw=weather_data[0]['WindSpeed [mph]']
        data_poly=wind_poly.iloc[:,0]
        set_value('coeff##cubic',wind_param.iloc[0,0:3].tolist())
        set_value('r2##cubic',wind_param.iloc[0,4].tolist())
    elif data=='Day 2':
        data_raw=weather_data[1]['WindSpeed [mph]']
        data_poly=wind_poly.iloc[:,1]
        set_value('coeff##cubic',wind_param.iloc[1,0:3].tolist())
        set_value('r2##cubic',wind_param.iloc[1,4].tolist())
    elif data=='Day 3':
        data_raw=weather_data[2]['WindSpeed [mph]']
        data_poly=wind_poly.iloc[:,2]
        set_value('coeff##cubic',wind_param.iloc[2,0:3].tolist())
        set_value('r2##cubic',wind_param.iloc[2,4].tolist())

    data_x=weather_data[0]['Time [h]']
    add_line_series('windVStime','Raw Data', data_x.tolist(), data_raw.tolist(),weight=2)
    add_scatter_series('windVStime','Cubic Fit', data_x.tolist(), data_poly.tolist())
    
    show_item('windVStime')
    show_item('cubic_param')
    show_item('cubic_r2')
    show_item('Next##2')

def hide2(sender, data):
    hide_item('Wind Speed Approximation')
    show_item('Weather Data Input')

def show2(sender, data):

    wind_poly=pd.read_excel(path+'/RESULTS/wind_poly.xlsx',header=0,index_col=0)

    days_to_compute=int(get_value('days_compute'))
    days=[*range(0,days_to_compute)]
    weather_data=pd.read_excel(path+'/DATA/weather.xlsx',sheet_name=days,header=0)

    temp_atm,temp_dewpt,solar_rad,wind_speed=conversion_sewing(days_to_compute,weather_data,wind_poly)

    sew_export=np.r_[[temp_atm],[temp_dewpt],[solar_rad],[wind_speed]]

    col_data=['Atmospheric Temp [ºC]', 'DewPoint Temp [ºC]', 'Solar Radiation [W/m²]','Wind Speed [km/h]']

    exportftemp=pd.DataFrame(sew_export.transpose(),columns=col_data)
    filepath = path+'/RESULTS/weather_sewing.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    hide_item('Wind Speed Approximation')
    show_item('Unit Conversion + Sewing')

#############################################################################

def hide3(sender, data):
    hide_item('Unit Conversion + Sewing')
    show_item('Wind Speed Approximation')

def weather_data_conv_sew(sender, data):

    clear_plot('weatherVStime')

    sew_import=pd.read_excel(path+'/RESULTS/weather_sewing.xlsx',header=0)

    time=sew_import.iloc[:,0]
    temp_atm=sew_import.iloc[:,1]
    temp_dewpt=sew_import.iloc[:,2]
    solar_rad=sew_import.iloc[:,3]
    wind_speed=sew_import.iloc[:,4]

    if data=='Atmospheric Temp':
        add_line_series('weatherVStime','Atmosperic Temp [ºC]',time.tolist(),temp_atm.tolist(),weight=2)
    if data=='DewPoint Temp':
        add_line_series('weatherVStime','DewPoint Temp [ºC]',time.tolist(),temp_dewpt.tolist(),weight=2)
    if data=='Solar Radiation':
        add_line_series('weatherVStime','Solar Radiation [W/m²]',time.tolist(),solar_rad.tolist(),weight=2)
    if data=='Wind Speed':
        add_line_series('weatherVStime','Wind Speed [km/h]',time.tolist(),wind_speed.tolist(),weight=2)

    show_item('weatherVStime')
    show_item('Next##3')

def show3(sender, data):
    hide_item('Unit Conversion + Sewing')
    show_item('Thermal Profile Calc. Param. 1/2')

#############################################################################

def hide3_2(sender, data):
    hide_item('Thermal Profile Calc. Param. 1/2')
    show_item('Unit Conversion + Sewing')

def show3_2(sender, data):

    time_step=get_value('delta_t') # [s]
    node_spacing=get_value('delta_x')

    sew_import=pd.read_excel(path+'/RESULTS/weather_sewing.xlsx',header=0)
    
    time=sew_import.iloc[:,0].to_numpy()
    time_sec=time*3600 # [s]
    time_interp= np.arange(0,time_sec[-1]+time_step,time_step) # [s]

    temp_atm=sew_import.iloc[:,1].to_numpy()
    temp_dewpt=sew_import.iloc[:,2].to_numpy()
    solar_rad=sew_import.iloc[:,3].to_numpy()
    wind_speed=sew_import.iloc[:,4].to_numpy()

    temp_atm_interp = np.interp(time_interp,time_sec,temp_atm)
    temp_dewpt_interp = np.interp(time_interp,time_sec,temp_dewpt)
    solar_rad_interp = np.interp(time_interp,time_sec,solar_rad)
    wind_speed_interp = np.interp(time_interp,time_sec,wind_speed)

    interp_export=np.c_[temp_atm_interp,temp_dewpt_interp,solar_rad_interp,wind_speed_interp]

    col_data=['Atmospheric Temp [ºC]', 'DewPoint Temp [ºC]', 'Solar Radiation [W/m²]','Wind Speed [km/h]']

    exportftemp=pd.DataFrame(interp_export,columns=col_data,index=time_interp)
    filepath = path+'/RESULTS/weather_interp.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    number_of_layers=int(get_value('number_layers'))

    if number_of_layers==2:
        show_item('text_1')
        hide_item('text_2')
        hide_item('text_3')
        hide_item('text_4')
        #
        show_item('density_1')
        hide_item('density_2')
        hide_item('density_3')
        hide_item('density_4')
        #
        show_item('specheat_1')
        hide_item('specheat_2')
        hide_item('specheat_3')
        hide_item('specheat_4')
        #
        show_item('conduc_1')
        hide_item('conduc_2')
        hide_item('conduc_3')
        hide_item('conduc_4')
        #
        show_item('thick_1')
        hide_item('thick_2')
        hide_item('thick_3')
        hide_item('thick_4')
        #
        show_item('contactR_1')
        hide_item('contactR_2')
        hide_item('contactR_3')
        hide_item('contactR_4')
    elif number_of_layers==3:
        hide_item('text_1')
        show_item('text_2')
        hide_item('text_3')
        hide_item('text_4')
        #
        hide_item('density_1')
        show_item('density_2')
        hide_item('density_3')
        hide_item('density_4')
        #
        hide_item('specheat_1')
        show_item('specheat_2')
        hide_item('specheat_3')
        hide_item('specheat_4')
        #
        hide_item('conduc_1')
        show_item('conduc_2')
        hide_item('conduc_3')
        hide_item('conduc_4')
        #
        hide_item('thick_1')
        show_item('thick_2')
        hide_item('thick_3')
        hide_item('thick_4')
        #
        hide_item('contactR_1')
        show_item('contactR_2')
        hide_item('contactR_3')
        hide_item('contactR_4')
    elif number_of_layers==4:
        hide_item('text_1')
        hide_item('text_2')
        show_item('text_3')
        hide_item('text_4')
        #
        hide_item('density_1')
        hide_item('density_2')
        show_item('density_3')
        hide_item('density_4')
        #
        hide_item('specheat_1')
        hide_item('specheat_2')
        show_item('specheat_3')
        hide_item('specheat_4')
        #
        hide_item('conduc_1')
        hide_item('conduc_2')
        show_item('conduc_3')
        hide_item('conduc_4')
        #
        hide_item('thick_1')
        hide_item('thick_2')
        show_item('thick_3')
        hide_item('thick_4')
        #
        hide_item('contactR_1')
        hide_item('contactR_2')
        show_item('contactR_3')
        hide_item('contactR_4')
    elif number_of_layers==5:
        hide_item('text_1')
        hide_item('text_2')
        hide_item('text_3')
        show_item('text_4')
        #
        hide_item('density_1')
        hide_item('density_2')
        hide_item('density_3')
        show_item('density_4')
        #
        hide_item('specheat_1')
        hide_item('specheat_2')
        hide_item('specheat_3')
        show_item('specheat_4')
        #
        hide_item('conduc_1')
        hide_item('conduc_2')
        hide_item('conduc_3')
        show_item('conduc_4')
        #
        hide_item('thick_1')
        hide_item('thick_2')
        hide_item('thick_3')
        show_item('thick_4')
        #
        hide_item('contactR_1')
        hide_item('contactR_2')
        hide_item('contactR_3')
        show_item('contactR_4')

    hide_item('Thermal Profile Calc. Param. 1/2')
    show_item('Thermal Profile Calc. Param. 2/2')

#############################################################################

def hide3_3(sender, data):
    hide_item('Thermal Profile Calc. Param. 2/2')
    show_item('Thermal Profile Calc. Param. 1/2')

def show3_3(sender, data):

    hide_item('Thermal Profile Calc. Param. 2/2')
    show_item('Thermal Profile Calculation Process')

    set_value('ite_number','')
    set_value('norm','')
    set_value('event','')
    set_value('event2','')
    set_value('event3','')
    hide_item('Next##34')

    number_of_layers = int(get_value('number_layers'))
    time_step = get_value('delta_t')

    node_spacing = get_value('delta_x')

    numite = get_value('numite')
    convcrit = get_value('convcrit')

    ground_density = get_value('density_ground')
    ground_specheat = get_value('specheat_ground')
    ground_conduc = get_value('conduc_ground')

    if number_of_layers==2:
        density = [get_value('density_1')]
        spec_heat_capacity = [get_value('specheat_1')]
        thermal_conductivity = [get_value('conduc_1')]
        pavement_thickness = [get_value('thick_1')]
        interface_contact_res = [get_value('contactR_1')]
    elif number_of_layers==3:
        density = get_value('density_2')
        spec_heat_capacity = get_value('specheat_2')
        thermal_conductivity = get_value('conduc_2')
        pavement_thickness = get_value('thick_2')
        interface_contact_res = get_value('contactR_2')
    elif number_of_layers==4:
        density = get_value('density_3')
        spec_heat_capacity = get_value('specheat_3')
        thermal_conductivity = get_value('conduc_3')
        pavement_thickness = get_value('thick_3')
        interface_contact_res = get_value('contactR_3')
    elif number_of_layers==4:
        density = get_value('density_4')
        spec_heat_capacity = get_value('specheat_4')
        thermal_conductivity = get_value('conduc_4')
        pavement_thickness = get_value('thick_4')
        interface_contact_res = get_value('contactR_4')

    pavement_thickness = np.round(pavement_thickness,4)

    ground_layer_depth = get_value('ground_layer')
    deep_ground_temp = get_value('ground_temp')+273.15

    surface_albedo = get_value('albedo')
    surface_emissivity = get_value('emissivity')
    sky_view_factor = get_value('sky_viewfactor')
    solar_view_factor = get_value('solar_viewfactor')
    characteristic_length = get_value('charac_length')

    simul_prop=[number_of_layers,time_step,node_spacing]
    pavement_depth = pavement_depth_const(pavement_thickness)
    ground_prop = [ground_layer_depth,deep_ground_temp]
    surface_prop=[surface_albedo,surface_emissivity,sky_view_factor,solar_view_factor,characteristic_length]

    density.append(ground_density)
    spec_heat_capacity.append(ground_specheat)
    thermal_conductivity.append(ground_conduc)

    interp_import = pd.read_excel(path+'/RESULTS/weather_interp.xlsx',header=0)

    profile,normM,max_ite,unstab_flag,convergence_flag = thermal_profile_calc(simul_prop,pavement_depth,ground_prop,surface_prop, \
                    density,spec_heat_capacity,thermal_conductivity,interface_contact_res, \
                    interp_import, \
                    stefan_boltzmann,ref_temp_flow,ref_viscosity_flow,ref_conductivity_flow,ref_diffusivity_flow,ref_prandtl_flow, \
                    numite, convcrit)

    normM_calc = normM[:max_ite+1]
    ite_order=np.arange(1,max_ite+2)

    depth_col = np.arange(0,(int(ground_layer_depth/node_spacing)+1)*node_spacing,node_spacing)

    '''
    exportftemp=pd.DataFrame(profile,index=interp_import.iloc[:,0],columns=depth_col)
    filepath = path+'/RESULTS/temp_profile.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)
    '''

    profile_celsius=profile-273.15
    
    profile_pd = pd.DataFrame(profile_celsius)
    profile_pd.to_csv(path+'/RESULTS/temp_profile.csv')

    profile_sub = profile_celsius[::10,:]

    exportftemp=pd.DataFrame(profile_sub,index=interp_import.iloc[::10,0],columns=depth_col)
    filepath = path+'/RESULTS/temp_profile_subsample.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    exportftemp=pd.DataFrame(normM_calc.transpose(),index=ite_order,columns=['Frobenius Norm'])
    filepath = path+'/RESULTS/recurring_iterations_log.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    if (convergence_flag == False) and (unstab_flag == False):
        set_value('event','ITERATION PROCESS ENDED W/O REACHING CONVERGENCE...')
        set_value('event2','Use the result as is (Next)')
        set_value('event3','OR change iteration Nº and/or convergence crit. (Recalc)')
        configure_item('numite',enabled=True)
        configure_item('convcrit',enabled=True)
        show_item('Recalculate')

    if unstab_flag == False:
        show_item('Next##34')

#############################################################################

def hide3_4(sender, data):
    hide_item('Thermal Profile Calculation Process')
    show_item('Thermal Profile Calc. Param. 2/2')

def show3_4(sender, data):
    hide_item('Thermal Profile Calculation Process')
    show_item('Thermal Profile / Time')

#############################################################################

def hide4(sender, data):
    hide_item('Thermal Profile / Time')
    show_item('Thermal Profile Calculation Process')
    clear_plot('tempVStime')

def plot_thprofile(sender, data):
    clear_plot('tempVStime')

    temp_profile = pd.read_excel(path+'/RESULTS/temp_profile_subsample.xlsx',header=1)
    time = temp_profile.iloc[:,0]/3600
    timemin = np.round(time*60,0)
    timemin = timemin.tolist()
    timeround = np.round(time,2)
    timeround = timeround.tolist()
    time = time.tolist()

    number_of_layers = int(get_value('number_layers'))
    node_spacing = get_value('delta_x')

    if number_of_layers==2:
        pavement_thickness = [get_value('thick_1')]
    elif number_of_layers==3:
        pavement_thickness = get_value('thick_2')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_3')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_4')

    pavement_thickness = np.round(pavement_thickness,4)
    pavement_depth = pavement_depth_const(pavement_thickness)
    pave = int(pavement_depth[-1]/node_spacing)

    for k in range(1,pave+2):
        temp = temp_profile.iloc[:,k]
        temp = temp.tolist()

        add_line_series('tempVStime', 'Depth = '+str(round((k-1)*node_spacing*100,2))+' cm', time, temp)

    configure_item('time_sel',items=timeround)
    configure_item('time_sel##2',items=timeround)

    show_item('Next##4')

def show4(sender, data):
    hide_item('Thermal Profile / Time')
    show_item('Thermal Profile / Depth')

#############################################################################

def hide5(sender, data):
    hide_item('Thermal Profile / Depth')
    show_item('Thermal Profile / Time')
    hide_item('Next##5')

    clear_plot('tempVSdepth')

def depthVStime(sender, data):
    clear_plot('tempVSdepth')

    temp_profile=pd.read_excel(path+'/RESULTS/temp_profile_subsample.xlsx',header=1)
    time=temp_profile.iloc[:,0]/3600
    time=np.round(time,2)

    node_spacing = get_value('delta_x')
    ground_layer_depth = get_value('ground_layer')

    vecDepth=np.arange(0,ground_layer_depth,node_spacing)

    tsel=time[time==float(get_value('time_sel'))].index[0]

    temp=temp_profile.iloc[tsel,1:]
    temp=temp.tolist()

    vecDepth_inv = - vecDepth
    vecDepth_inv = vecDepth_inv.tolist()

    add_line_series('tempVSdepth','', temp, vecDepth_inv, weight=2, color=[255, 0, 0, 255])
    configure_item('tempVSdepth',label='Temperature Profile - '+get_value('time_sel')+' h')

    show_item('Export')

    calc_choice = get_value('calc_choice')

    if calc_choice == 'Thermal + Stress Profile' :
        show_item('Next##5')

def export(sender,data):

    temp_profile=pd.read_excel(path+'/RESULTS/temp_profile_subsample.xlsx',header=1)
    time=temp_profile.iloc[:,0]/3600
    time=np.round(time,2)

    node_spacing = get_value('delta_x')
    ground_layer_depth = get_value('ground_layer')

    vecDepth=np.arange(0,ground_layer_depth,node_spacing)

    tsel=time[time==float(get_value('time_sel'))].index[0]
    tsel_2=time[time==float(get_value('time_sel'))].tolist()

    exportftemp=pd.DataFrame(temp_profile.iloc[tsel,1:].tolist(),index=[vecDepth],columns=[tsel_2[0]])
    filepath = path+'/RESULTS/depthVStemp_'+str(tsel_2[0])+'h.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

def show5(sender, data):
    hide_item('Thermal Profile / Depth')
    show_item('Thermal Contraction Coeff')

#############################################################################

def hide5_1(sender, data):
    hide_item('Thermal Contraction Coeff')
    show_item('Thermal Profile / Depth')

def coeff_sel(sender, data):
    if get_value('thermal_coeff')==combo_item[0]:
        show_item('direct_in')
        hide_item('partial_in')
        hide_item('separate_in')
    elif get_value('thermal_coeff')==combo_item[1]:
        hide_item('direct_in')
        show_item('partial_in')
        hide_item('separate_in')
    elif get_value('thermal_coeff')==combo_item[2]:
        hide_item('direct_in')
        hide_item('partial_in')
        show_item('separate_in')
    else:
        hide_item('direct_in')
        hide_item('partial_in')
        hide_item('separate_in')
    hide_item('Next##51')

def valid_coeff(sender, data):
    if get_value('thermal_coeff')==combo_item[0]:
        bmix=get_value('bmix')
    elif get_value('thermal_coeff')==combo_item[1]:
        lbin=get_value('lbin##1')
        lag=get_value('lag##1')
        vma=get_value('vma')
        veff=get_value('veff')
        bmix=(vma*lbin+veff*lag)/300
    elif get_value('thermal_coeff')==combo_item[2]:
        lbin=get_value('lbin##2')
        lag=get_value('lag##2')
        gsb=get_value('gsb')
        gb=get_value('gb')
        gmm=get_value('gmm')
        av=get_value('av')
        bc=get_value('bc')
        gmb=(1-av/100)*gmm
        vma=round(100-gmb*(100-bc)/gsb,2)
        gse=(100-bc)/(100/gmm-bc/gb)
        veff=round(gmb*(100-bc)/gse,2)
        bmix=(vma*lbin+veff*lag)/300
    set_value('bmix_value', bmix)
    show_item('bmix_val')
    show_item('Next##51')

def show5_1(sender, data):
    hide_item('Thermal Contraction Coeff')
    show_item('Mechanical Properties Input')

#############################################################################

def hide6(sender, data):
    hide_item('Mechanical Properties Input')
    show_item('Thermal Profile / Depth')
    clear_table('mat_table')
    hide_item('Next##6')

def mat_input(sender,data):

    datain=get_value('mat_choice')

    if datain==mat_choice[0]:
        col_name=['Temp [ºC]','Time [s]','Compliance [1/GPa]']
        datain=pd.read_excel(path+'/DATA/creep_compliance.xlsx')
        configure_item('Experimental Data',x_axis_name='Time [s]', \
                       y_axis_name='Creep Compliance [1/GPa]')
        configure_item('Master Curve##plot',x_axis_name='Reduced Time [s]', \
                       y_axis_name='Creep Compliance [1/GPa]')

    elif datain==mat_choice[1]:
        # col_name=['Temp [ºC]','Time [s]','Erelax [GPa]']
        col_name=['Temp [ºC]','Freq [Hz]','Erelax [GPa]']
        datain=pd.read_excel(path+'/DATA/relax_modulus.xlsx')
        configure_item('Experimental Data',x_axis_name='Time [s]', \
                       y_axis_name='Relaxation Modulus [GPa]')
        configure_item('Master Curve##plot',x_axis_name='Reduced Time [s]', \
                       y_axis_name='Relaxation Modulus [GPa]')

    table_in=datain.values
    table_in=table_in.tolist()

    clear_table('mat_table')
    set_headers('mat_table',col_name)
    set_table_data('mat_table',table_in)
    show_item('mat_table')
    show_item('Next##6')

def show6(sender, data):
    hide_item('Mechanical Properties Input')
    show_item('Log-Log Optim Experimental Data')

#############################################################################

def hide7(sender, data):
    hide_item('Log-Log Optim Experimental Data')
    show_item('Mechanical Properties Input')
    clear_plot('Experimental Data')
    hide_item('Next##7')
    hide_item('Optimization Poly(log-log)')
    hide_item('Export##2')

def plot_expe(sender,data):
    clear_plot('Experimental Data')

    temps = [float(row[0]) for row in data]
    datain = get_value('mat_choice')
    time = [float(row[1]) for row in data]
    if datain==mat_choice[1]:
        time = np.array(time)
        time = 1/time
        time = time.tolist()
    val = [float(row[2]) for row in data]

    temp_list=list(temps)
    temp_set=set(temps)
    temp_set=sorted(temp_set,key=temp_list.index)
    temp_set=list(temp_set)

    if sender=='Plot Data':
        configure_item('Experimental Data',xaxis_log_scale=False,yaxis_log_scale=False)
    elif sender=='Plot Log-Log Data':
        configure_item('Experimental Data',xaxis_log_scale=True,yaxis_log_scale=True)

    for tu in range(len(temp_set)):
        data_x=[]
        data_y=[]
        for k in range(len(temps)):
            if temps[k]==temp_set[tu]:
                data_x.append(time[k])
                data_y.append(val[k])
        add_line_series('Experimental Data',str(temp_set[tu])+' ºC', data_x, data_y)

    hide_item('Next##7')
    show_item('Optimization Poly(log-log)')
    configure_item('tref##sel',items=temp_set)

def optim_expe(sender,data):

    material_data=pd.DataFrame(data,dtype=float)
    temperature_data=material_data.iloc[:,0].unique()

    if get_value('mat_choice')==mat_choice[0]:
        typeI=1
        configure_item('Experimental Data',x_axis_name='Time [s]',y_axis_name='Creep Compliance [1/GPa]')
    else:
        typeI=0
        configure_item('Experimental Data',x_axis_name='Time [s]',y_axis_name='Relaxation Modulus [GPa]')

    raw_coeff,raw_stderr=optim_loglog(material_data,typeI,temperature_data)

    configure_item('Experimental Data',xaxis_log_scale=True,yaxis_log_scale=True)
    configure_item('Experimental Data',label='Optimized Curves')

    i=0

    if sender=='Export##2': text_file = open(path+'/RESULTS/polymial_optim_material.txt', 'w')
    
    for k in temperature_data:

        tempselection=material_data[material_data.iloc[:,0]==k]
        mat=np.log10(tempselection.iloc[:,2])

        if typeI==1:
            t=np.log10(tempselection.iloc[:,1])
            mat_calc=polylogJ(t,*raw_coeff[i,:])
        else:
            t=np.log10(1/tempselection.iloc[:,1])
            mat_calc=linlogE(t,*raw_coeff[i,:])
        
        if sender=='Optimization Poly(log-log)':
            data_x=[]
            data_y=[]
            for n in range(len(mat_calc)):
                data_x.append(10**(t.iloc[n]))
                data_y.append(10**(mat_calc.iloc[n]))
            add_scatter_series('Experimental Data',str(temperature_data[i])+' ºC', data_x, data_y)
            show_item('Next##7')
            show_item('Export##2')

        elif sender=='Export##2':
            if typeI==1:
                text_file.write('Polynomial fit "a.x²+b.x+c" for temperature = '+str(k)+' ºC\nResults :')
                text_file.write('\na = '+str(raw_coeff[i,0])+' +/- '+str(raw_stderr[i,0,0]**0.5))
                text_file.write('\nb = '+str(raw_coeff[i,1])+' +/- '+str(raw_stderr[i,1,1]**0.5))
                text_file.write('\nc = '+str(raw_coeff[i,2])+' +/- '+str(raw_stderr[i,2,2]**0.5))
            else:
                text_file.write('Linear fit "a.x+b" for temperature ='+str(k)+'ºC\nResults :')
                text_file.write('\na = '+str(raw_coeff[i,0])+' +/- '+str(raw_stderr[i,0,0]**0.5))
                text_file.write('\nb = '+str(raw_coeff[i,1])+' +/- '+str(raw_stderr[i,1,1]**0.5))

            residuals=mat-mat_calc
            ss_res=np.sum(residuals**2)
            ss_tot=np.sum((mat-np.mean(mat))**2)
            r_squared=1-(ss_res/ss_tot)
            text_file.write('\nR² = '+str(r_squared))
            text_file.write('\n')
            text_file.write('\n')

        i+=1

    if sender=='Export##2': text_file.close()

def show7(sender, data):
    hide_item('Log-Log Optim Experimental Data')
    show_item('Master Curve')

#############################################################################
    
def hide8(sender, data):
    hide_item('Master Curve')
    show_item('Log-Log Optim Experimental Data')
    hide_item('Next##8')
    clear_plot('Master Curve##plot')

def tref(sender,data):
    configure_item('Build Master Curve',enabled=True)
    hide_item('Next##8')

def master_curve(sender,data):
    temps = [float(row[0]) for row in data]
    temp_list=list(temps)
    temp_set=set(temps)
    temp_set=sorted(temp_set,key=temp_list.index)
    temp_set=np.array(temp_set)

    tref_chosen=float(get_value('tref##sel'))

    trefInd=np.where(temp_set==tref_chosen)
    tref_index=int(trefInd[0])
    temp_minus_tref=np.delete(temp_set,trefInd)

    material_data=pd.DataFrame(data,dtype=float)
    temperature_data=material_data.iloc[:,0].unique()

    if get_value('mat_choice')==mat_choice[0]:
        typeI=1
        name='creep_compliance'
        configure_item('shift_factor',items=['William-Landel-Ferry','Arrhenius'])
    else:
        typeI=0
        name='relax_modulus'
        configure_item('shift_factor',items=['Arrhenius'])

    raw_coeff,raw_stderr=optim_loglog(material_data,typeI,temperature_data)
    
    shifting_coeff=mean_shift(material_data,typeI,temp_minus_tref,tref_index,raw_coeff)

    reduced_time=equiv_slope(material_data,typeI,temp_minus_tref,tref_chosen,shifting_coeff)

    tr_order,MC_order=reduced_t_order(reduced_time,material_data)

    data_x=[]
    data_y=[]

    for i in range(len(tr_order)):
        data_x.append(tr_order[i])
        data_y.append(MC_order[i])
    
    add_scatter_series('Master Curve##plot','', data_x, data_y)

    coeff_wlf,stderr_wlf=optim_wlf(shifting_coeff,temperature_data,tref_chosen,tref_index)
    coeff_arrh,stderr_arrh=optim_arrhenius(shifting_coeff,temperature_data,tref_chosen,tref_index)

    shifting_coeff=np.insert(shifting_coeff,tref_index,0)

    exportftemp=pd.DataFrame(shifting_coeff,index=[temp_set])
    filepath = path+'/RESULTS/discrete-shift.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    exportftemp=pd.DataFrame(MC_order,index=[tr_order])
    filepath = path+'/RESULTS/master-curve_'+name+'.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    exportftemp=pd.DataFrame([[coeff_wlf[0],stderr_wlf[0,0]**0.5], \
                              [coeff_wlf[1],stderr_wlf[1,1]**0.5]], \
                              columns=['coeff','stderr'],index=['C1','C2'])
    filepath = path+'/RESULTS/optim-wlf.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    exportftemp=pd.DataFrame([[coeff_arrh[0],stderr_arrh[0,0]**0.5]], \
                              columns=['coeff','stderr'],index=['Ea'])
    filepath = path+'/RESULTS/optim-arrhenius.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    show_item('Next##8')
    
def show8(sender, data):
    hide_item('Master Curve')
    show_item('Shift-Factor')

#############################################################################

def hide9(sender,data):
    hide_item('Shift-Factor')
    show_item('Master Curve')
    configure_item('shift_plot',label='')
    hide_item('wlf')
    hide_item('arrh')
    hide_item('r2')
    hide_item('Next##9')
    clear_plot('shift_plot')

def shift(sender,data):
    shift_chosen=get_value('shift_factor')
    shifting_coeff=pd.read_excel(path+'/RESULTS/discrete-shift.xlsx',header=0)

    shifting_coeff=shifting_coeff.iloc[:,1].to_numpy()

    tref_chosen=float(get_value('tref##sel'))

    if shift_chosen=='William-Landel-Ferry':
        coeff=pd.read_excel(path+'/RESULTS/optim-wlf.xlsx',header=0)
        fit_law=make_wlf(tref_chosen)
        configure_item('shift_plot',label='William-Landel-Ferry')
        coeff=coeff.iloc[:,1].to_numpy()
        set_value('coeff##wlf',coeff)
        show_item('wlf')
        hide_item('arrh')
    elif shift_chosen=='Arrhenius':
        coeff=pd.read_excel(path+'/RESULTS/optim-arrhenius.xlsx',header=0)
        configure_item('shift_plot',label='Arrhenius')
        fit_law=make_arrhenius(tref_chosen)
        coeff=coeff.iloc[:,1].to_numpy()
        set_value('coeff##arrh',coeff[0])
        show_item('arrh')
        hide_item('wlf')

    temps = [float(row[0]) for row in data]
    temp_list=list(temps)
    temp_set=set(temps)
    temp_set=sorted(temp_set,key=temp_list.index)
    temp_set=np.array(temp_set)

    logat_calc=fit_law(temp_set,*coeff)

    residuals = shifting_coeff- logat_calc
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((shifting_coeff-np.mean(shifting_coeff))**2)
    r_squared = 1 - (ss_res / ss_tot)
    set_value('r_squared',r_squared)

    temprange=np.arange(-50,80,1)
    logat_sim=fit_law(temprange,*coeff)

    add_scatter_series('shift_plot','Discrete Values',temp_set,shifting_coeff)
    add_line_series('shift_plot','Optimized Law',temprange,logat_sim)

    show_item('r2')
    show_item('Next##9')

def show9(sender,data):
    hide_item('Shift-Factor')
    if get_value('mat_choice')==mat_choice[0]:
        show_item('Prony Optimization##1')
    else:
        show_item('Prony Optimization##2')


#############################################################################

def hide10(sender,data):
    hide_item('Prony Optimization##1')
    show_item('Shift-Factor')
    configure_item('CCMC Optim',enabled=False)
    hide_item('Interconversion')
    clear_plot('Material Optimization##plot')
    hide_item('Next##10')

def ccmc_optim(sender,data):
    
    if sender=='nn##CCMC':
        configure_item('CCMC Optim',enabled=True)
        hide_item('Interconversion')
        clear_plot('Material Optimization##plot')
        configure_item('Material Optimization##plot',yaxis2=False)
        configure_item('Material Optimization##plot',y_axis_name='Creep Compliance [1/GPa]')

    elif sender=='CCMC Optim':
        clear_plot('Material Optimization##plot')
        configure_item('Material Optimization##plot',yaxis2=False, \
                       y_axis_name='Creep Compliance [1/GPa]')

        master_data=pd.read_excel(path+'/RESULTS/master-curve_creep_compliance.xlsx',header=0)
        tr_order=master_data.iloc[:,0].to_numpy()
        MC_order=master_data.iloc[:,1].to_numpy()
        nn=int(get_value('nn##CCMC'))

        fit_gpower=make_GpowerLaw(nn)
        coeff,stderr=curve_fit(fit_gpower,tr_order,MC_order,p0=[1.0]*(2*nn+2),bounds=(0.0,np.inf),method='trf')

        index_list=[]
        index_list.append('D_0')
        for i in range(1,nn+1):
            index_list.append('D_'+str(i))
            index_list.append('tau_'+str(i))
        index_list.append('k')

        exportftemp=pd.DataFrame(coeff,index=[index_list])
        filepath = path+'/RESULTS/optim-GenModPower_CCMC.xlsx'
        exportftemp.to_excel(filepath,index=True,header=True)

        tsim=np.logspace(-6,12,50)
        mat_sim=fit_gpower(tsim,*coeff)

        add_scatter_series('Material Optimization##plot','CCMC Discrete Values',tr_order,MC_order)
        add_line_series('Material Optimization##plot','CCMC Optimized Law', \
                        tsim.tolist(),mat_sim.tolist(),weight=2.0)

        show_item('Interconversion')

    elif sender=='Interconversion':

        configure_item('Material Optimization##plot',yaxis2=True, \
                       y_axis_name='CCMC [1/GPa] (left) / ErelMC [GPa] (right)', \
                       y2axis_log_scale=True)
        
        coeff=pd.read_excel(path+'/RESULTS/optim-GenModPower_CCMC.xlsx',header=0)
        coeff=coeff.iloc[:,1].to_numpy()
        nn=int(get_value('nn##CCMC'))

        master_data=pd.read_excel(path+'/RESULTS/master-curve_creep_compliance.xlsx',header=0)
        tr=master_data.iloc[:,0].to_numpy()
        CCMC=master_data.iloc[:,1].to_numpy()
                
        fit_gpower=make_GpowerLaw(nn)
        tr_order=np.logspace(-6,12,50)
        MC_order=fit_gpower(tr_order,*coeff)
        
        Edyn_IC=general_interconv(tr_order,MC_order,coeff,nn,'GMpower')

        tr_orderE=tr_order[0:len(tr_order)-1]

        exportftemp=pd.DataFrame(Edyn_IC,index=[tr_orderE])
        filepath = path+'/RESULTS/interconversion_relax_modulus.xlsx'
        exportftemp.to_excel(filepath,index=True,header=True)

        add_line_series('Material Optimization##plot','Interc. Erel (right)', \
                        tr_orderE,Edyn_IC,axis=1,weight=2.0)

        show_item('Next##10')

def show10(sender,data):
    hide_item('Prony Optimization##1')
    show_item('Prony Optimization##2')

#############################################################################

def hide11(sender,data):
    hide_item('Prony Optimization##2')
    clear_plot('Material Optimization##plot2')
    configure_item('ErelMC Optim',enabled=False)
    hide_item('Next##11')
    if get_value('mat_choice')==mat_choice[0]:
        show_item('Prony Optimization##1')
        configure_item('CCMC Optim',enabled=False)
        hide_item('Interconversion')
        clear_plot('Material Optimization##plot')
        hide_item('Next##10')
    else:
        show_item('Shift-Factor')

def prony_optim2(sender,data):
    flag_error = False

    if sender=='nn##ErelMC':
        configure_item('ErelMC Optim',enabled=True)
        clear_plot('Material Optimization##plot2')

    elif sender=='ErelMC Optim':
        nn=int(get_value('nn##ErelMC'))

        if get_value('mat_choice')==mat_choice[0]:
            master_data=pd.read_excel(path+'/RESULTS/interconversion_relax_modulus.xlsx',header=0)
            tr_order=master_data.iloc[:,0].to_numpy()
            MC_order=master_data.iloc[:,1].to_numpy()
            coeff_prony,stderr_prony=optim_interconv(tr_order,MC_order,nn)
        else:
            master_data=pd.read_excel(path+'/RESULTS/master-curve_relax_modulus.xlsx',header=0)
            tr_order=master_data.iloc[:,0].to_numpy()
            MC_order=master_data.iloc[:,1].to_numpy()
            typeI=0
            try:
                #coeff_sigm,stderr_sigm=curve_fit(sigmoid,tr_order,np.log10(MC_order),p0=[0.0]*4,method='lm') # methods : lm / trf / dogbox
                x0 = np.array([0,0,0,0])
                res = least_squares(sigmoidLSQ,x0,args=(tr_order,np.log10(MC_order)))
                #res = minimize(sigmoidLSQ,x0,args=(tr_order,np.log10(MC_order)),method='trust-constr')
                coeff_sigm = res.x
                trsim=np.logspace(-6,12,50)
                logErelSim=sigmoid(trsim,coeff_sigm)
                ErelSim=10**logErelSim
                coeff_prony,stderr_prony=optim_prony(trsim,ErelSim,typeI,nn)
            except:
                coeff_prony,stderr_prony=optim_prony(tr_order,MC_order,typeI,nn)
                flag_error = True

        index_list=[]
        index_list.append('E_0')
        for i in range(1,nn+1):
            index_list.append('E_'+str(i))
            index_list.append('log(tau_'+str(i)+')')

        exportftemp=pd.DataFrame(coeff_prony,index=[index_list])
        filepath = path+'/RESULTS/optim-prony_ErelMC.xlsx'
        exportftemp.to_excel(filepath,index=True,header=True)

        fit_erelmc=make_pronyE(nn)
        tsim=np.logspace(-6,12,50)
        mat_sim=fit_erelmc(tsim,*coeff_prony)

        add_scatter_series('Material Optimization##plot2','ErelMC Discrete Values',tr_order,MC_order)
        
        if  get_value('mat_choice')==mat_choice[1] and not(flag_error):
            add_line_series('Material Optimization##plot2','ErelMC Optimized Sigmoid', \
                            trsim.tolist(),ErelSim.tolist())
        

        add_line_series('Material Optimization##plot2','ErelMC Optimized Prony', \
                        tsim.tolist(),mat_sim.tolist(),weight=2.0)

        show_item('Next##11')

def show11(sender,data):
    hide_item('Prony Optimization##2')

    temp_profile=pd.read_excel(path+'/RESULTS/temp_profile_subsample.xlsx',header=1)
    #temp_profile=pd.read_excel(path+'/RESULTS/temp_profile.xlsx',header=1)

    vecTime=temp_profile.iloc[:,0].to_numpy()

    number_of_layers = int(get_value('number_layers'))
    node_spacing = get_value('delta_x')
    ground_layer_depth = get_value('ground_layer')

    if number_of_layers==2:
        pavement_thickness = [get_value('thick_1')]
    elif number_of_layers==3:
        pavement_thickness = get_value('thick_2')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_3')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_4')

    pavement_thickness = np.round(pavement_thickness,4)
    pave = int(pavement_thickness[0]/node_spacing)

    vecDepth=np.arange(0,ground_layer_depth,node_spacing)
    vecDepth = vecDepth[0:pave+1]

    temp_profile=temp_profile.iloc[:,1:(pave+2)].to_numpy()
    temp_profile=temp_profile+273.15

    #vecTimeInterp=timeinterp(vecTime)
    #vecTempInterp=tempinterp(vecTime,vecTimeInterp,temp_profile,vecDepth)
    vecTimeInterp=vecTime
    vecTempInterp=temp_profile

    shift_chosen=get_value('shift_factor')
    tref=float(get_value('tref##sel'))
    bmix=float(get_value('bmix_value'))
    #tmean=float(get_value('tmean'))
    tmean=get_value('ground_temp')+273.15 # equilibrium temperature

    if shift_chosen=='William-Landel-Ferry':
        coeff=pd.read_excel(path+'/RESULTS/optim-wlf.xlsx',header=0)
        coeff=coeff.iloc[:,1].to_numpy()
        logat=tsecWLF(vecTempInterp,coeff,tref)
    elif shift_chosen=='Arrhenius':
        coeff=pd.read_excel(path+'/RESULTS/optim-arrhenius.xlsx',header=0)
        coeff=coeff.iloc[:,1].to_numpy()
        logat=tsecArrh(vecTempInterp,coeff,tref)

    RedTime,Delta_RedTime,Etot,Delta_Etot=stressprecalc(vecTimeInterp, \
                                          vecTempInterp,vecDepth,logat,bmix,tmean)

    Eref=1
    EpronyBranch=int(get_value('nn##ErelMC'))
    EpronyCoeff=pd.read_excel(path+'/RESULTS/optim-prony_ErelMC.xlsx',header=0)
    EpronyCoeff=EpronyCoeff.iloc[:,1].to_numpy()

    CalculatedStressPSEUDO=stresscalc(vecTimeInterp,vecDepth, \
                           Delta_RedTime,Etot,Delta_Etot, \
                           EpronyCoeff,EpronyBranch,Eref)

    '''
    CalculatedStressFD=stresscalcFD(vecTimeInterp,vecDepth, \
                       Delta_RedTime,Etot,Delta_Etot, \
                       EpronyCoeff,EpronyBranch)
    '''

    exportftemp=pd.DataFrame(CalculatedStressPSEUDO,index=[vecTime],columns=[vecDepth*100])
    filepath = path+'/RESULTS/calculated_stress_GPa.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

    show_item('Stress Profile / Time')

#############################################################################

def hide12(sender, data):
    hide_item('Stress Profile / Time')
    show_item('Prony Optimization##2')
    clear_plot('stressVStime')

def plot_stressprofile(sender, data):
    clear_plot('stressVStime')

    stress_profile=pd.read_excel(path+'/RESULTS/calculated_stress_GPa.xlsx',header=1)
    time=stress_profile.iloc[:,0]/3600
    time = time.tolist()

    number_of_layers = int(get_value('number_layers'))
    node_spacing = get_value('delta_x')

    if number_of_layers==2:
        pavement_thickness = [get_value('thick_1')]
    elif number_of_layers==3:
        pavement_thickness = get_value('thick_2')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_3')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_4')

    pavement_thickness = np.round(pavement_thickness,4)
    pave = int(pavement_thickness[0]/node_spacing)

    for k in range(1,pave+2):
        stress = stress_profile.iloc[:,k]*1e6         # [GPa] to [kPa]
        #stress = stress_profile.iloc[:,k]*145037.7    # [GPa] to [PSI]
        stress = stress.tolist()

        add_line_series('stressVStime','Depth = '+str(round((k-1)*node_spacing*100,2))+' cm', time, stress)

    show_item('Next##12')

def show12(sender, data):
    hide_item('Stress Profile / Time')
    show_item('Stress Profile / Depth')

#############################################################################

def hide13(sender, data):
    hide_item('Stress Profile / Depth')
    show_item('Stress Profile / Time')
    clear_plot('stressVSdepth')

def depthVSstress(sender, data):
    clear_plot('stressVSdepth')

    stress_profile=pd.read_excel(path+'/RESULTS/calculated_stress_GPa.xlsx',header=1)
    time=stress_profile.iloc[:,0]/3600
    time=np.round(time,2)

    number_of_layers = int(get_value('number_layers'))
    node_spacing = get_value('delta_x')
    ground_layer_depth = get_value('ground_layer')

    if number_of_layers==2:
        pavement_thickness = [get_value('thick_1')]
    elif number_of_layers==3:
        pavement_thickness = get_value('thick_2')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_3')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_4')

    pavement_thickness = np.round(pavement_thickness,4)
    pave = int(pavement_thickness[0]/node_spacing)

    vecDepth=np.arange(0,ground_layer_depth,node_spacing)
    vecDepth = vecDepth[0:pave+1]

    vecDepth_inv = - vecDepth
    vecDepth_inv = vecDepth_inv.tolist()

    tsel=time[time==float(get_value('time_sel##2'))].index[0]

    stress = stress_profile.iloc[tsel,1:]*1e6         # [GPa] to [kPa]
    #stress = stress_profile.iloc[tsel,1:]*145037.7    # [GPa] to [PSI]
    stress = stress.tolist()

    add_line_series('stressVSdepth','', stress, vecDepth_inv, weight=2, color=[0, 0, 0, 255])
    configure_item('stressVSdepth',label='Stress Profile of 1st Layer - '+get_value('time_sel##2')+' h')

    show_item('Export##3')

def exportstress(sender,data):

    stress_profile=pd.read_excel(path+'/RESULTS/calculated_stress_GPa.xlsx',header=1)
    time=stress_profile.iloc[:,0]/3600
    time=np.round(time,2)

    tsel=time[time==float(get_value('time_sel##2'))].index[0]
    tsel_2=time[time==float(get_value('time_sel##2'))].tolist()

    number_of_layers = int(get_value('number_layers'))
    node_spacing = get_value('delta_x')
    ground_layer_depth = get_value('ground_layer')

    if number_of_layers==2:
        pavement_thickness = [get_value('thick_1')]
    elif number_of_layers==3:
        pavement_thickness = get_value('thick_2')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_3')
    elif number_of_layers==4:
        pavement_thickness = get_value('thick_4')

    pavement_thickness = np.round(pavement_thickness,4)
    pave = int(pavement_thickness[0]/node_spacing)

    vecDepth=np.arange(0,ground_layer_depth,node_spacing)
    vecDepth = vecDepth[0:pave+1]

    exportftemp=pd.DataFrame(stress_profile.iloc[tsel,1:].tolist(),index=[vecDepth],columns=[tsel_2[0]])
    filepath = path+'/RESULTS/depthVSstress_GPa_'+str(tsel_2[0])+'h.xlsx'
    exportftemp.to_excel(filepath,index=True,header=True)

#############################################################################

with window('ACTS',no_resize=True,autosize=False):
   
    with tab_bar('tab'):

        with tab('Start Page'):
            add_spacing(count=2)
            add_same_line(spacing=20)
            add_button('ReadMe',callback=readme)
            add_same_line(spacing=440)
            add_button('About##button',callback=about)

            add_spacing(count=1)
            add_same_line(spacing=80)
            add_image(name='Logo',value=path+'/AC_THERMAL_STRESS/INTER_DATA/LOGO_NAME2.png')

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_text('Select what to compute:')
            add_combo('calc_choice',items=['Thermal + Stress Profile','Thermal Profile ONLY'], \
                      default_value='Thermal + Stress Profile', \
                      label='', width=200)

            add_spacing(count=1)
            add_same_line(spacing=510)
            add_button('START',show=True,callback=start)

        with tab('About',show=False):
            add_spacing(count=5)

            add_text('About the software:')
            add_separator()
            add_text('blablabla')

            add_spacing(count=5)
            add_same_line(spacing=280)
            add_button('BACK',show=True,callback=start2)

        with tab('Weather Data Input',show=False):
            add_spacing(count=5)

            add_combo('days_compute',items=['1','2','3'], \
                      label='Number of days to compute', width=50, \
                      callback=input_weather_data)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_combo('day_data##1',width=550,label='', \
                       callback=weather_test, \
                       callback_data=lambda:get_value('day_data##1'), \
                       show=False)
            
            add_spacing(count=2)

            add_table('weather_table',[],height=-45,show=False)

            add_spacing(count=5)
            add_same_line(spacing=20)
            add_button('StartPage',callback=start2)
            add_same_line(spacing=420)
            add_button('Next##1',show=False,callback=show)

        with tab('Wind Speed Approximation',show=False):
            add_spacing(count=2)

            add_combo('day_data##2',width=550,label='', \
                       callback=wind_speed_plot, \
                       callback_data=lambda:get_value('day_data##2'))

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_plot('windVStime',label='Wind Speed Cubic Polynomial Fit',height=270, \
                     x_axis_name='Time [h]',y_axis_name='Wind Speed [mph]',show=False)

            add_spacing(count=2)

            add_same_line(spacing=30)
            add_drag_float4(name='cubic_param',source='coeff##cubic',show=False, \
                            label=' [a,b,c,d] for ax³+bx²+cx+d', \
                            format='%.4f',width=300,no_input=True)
            
            add_spacing(count=2)

            add_same_line(spacing=258)
            add_drag_float(name='cubic_r2',source='r2##cubic',show=False, \
                           label=' = R²', \
                           format='%.4f',width=72,no_input=True)

            add_spacing(count=2)

            add_same_line(spacing=20)
            add_button('Previous##2',callback=hide2)
            add_same_line(spacing=440)
            add_button('Next##2',callback=show2,show=False)

        with tab('Unit Conversion + Sewing',show=False):
            add_spacing(count=2)

            add_combo('weather_choice',label='Sewed weather data', \
                      items=['Atmospheric Temp', 'DewPoint Temp', 'Solar Radiation','Wind Speed'], \
                      width=350,callback=weather_data_conv_sew, \
                      callback_data=lambda:get_value('weather_choice'))

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_plot('weatherVStime',label='Weather Data (all days)',height=315, \
                     x_axis_name='Time [h]',show=False)

            add_spacing(count=5)

            add_same_line(spacing=20)
            add_button('Previous##3',callback=hide3)
            add_same_line(spacing=440)
            add_button('Next##3',callback=show3,show=False)

        with tab('Thermal Profile Calc. Param. 1/2',show=False):
            add_spacing(count=2)

            add_drag_int(name='delta_t',label='Discrete Time-Step [s]',default_value=90,width=80)
            add_spacing(count=1)
            add_drag_float(name='delta_x',label='Discrete Spatial-Step [m] -- 0.0127 [m] ~ 0.5 [in]', \
                           default_value=0.0127,format='%.4f',width=80)
            
            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_text('Min. nº layers = 2 --> HMA + GROUND | Max. = 5 --> HMA + Old-HMA + ... + GROUND')
            add_spacing(count=2)
            add_combo(name='number_layers',label='Number of Layers of the slab [2-5]', \
                      items=['2','3','4','5'],default_value='2',width=80)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_text('SURFACE-MATERIAL PROPERTIES & CHARACTERISTICS')
            add_spacing(count=2)
            add_drag_float(name='albedo',label='Surface Albedo [0-1]',default_value=0.15,format='%.2f',width=80)
            add_drag_float(name='emissivity',label='Surface Infrared (IR) Emissivity [0-1]',default_value=0.97,format='%.2f',width=80)
            add_drag_float(name='sky_viewfactor',label='Sky View Factor [0-1]',default_value=0.95,format='%.2f',width=80)
            add_drag_float(name='solar_viewfactor',label='Solar View Factor [0-1]',default_value=0.85,format='%.2f',width=80)
            add_drag_float(name='charac_length',label='Characteristic Length / convection calc. [m]',default_value=1,format='%.2f',width=80)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_text('DEEP-GROUND PROPERTIES')
            add_spacing(count=2)
            add_drag_float(name='ground_temp',label='Deep Ground Temperature [ºC]',default_value=22,format='%.2f',width=80)
            add_drag_int(name='ground_layer',label='Max. Ground Depth [m]',default_value=3,width=80)

            add_spacing(count=3)

            add_same_line(spacing=20)
            add_button('Previous##32',callback=hide3_2)
            add_same_line(spacing=440)
            add_button('Next##32',callback=show3_2,show=True)

        with tab('Thermal Profile Calc. Param. 2/2',show=False):
            add_spacing(count=5)

            add_text('Material Properties for LAYER ... ')
            add_spacing(count=2)
            add_drag_int(name='text_1',label='',default_value=1,show=False,width=70,enabled=False)
            add_drag_int2(name='text_2',label='',default_value=[1,2],show=False,width=140,enabled=False)
            add_drag_int3(name='text_3',label='',default_value=[1,2,3],show=False,width=210,enabled=False)
            add_drag_int4(name='text_4',label='',default_value=[1,2,3,4],show=False,width=280,enabled=False)
            add_same_line(spacing=0)
            add_input_text(name='text_ground',label='',default_value='  GROUND',show=True,width=70,enabled=False)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_drag_float(name='density_1',default_value=2085,show=False,label='',format='%.0f',width=70)
            add_drag_float2(name='density_2',show=False,label='',format='%.0f',width=140)
            add_drag_float3(name='density_3',show=False,label='',format='%.0f',width=210)
            add_drag_float4(name='density_4',show=False,label='',format='%.0f',width=280)
            add_same_line(spacing=0)
            add_drag_float(name='density_ground',default_value=2200,format='%.0f',show=True,label='Density [kg/m³]',width=70)

            add_spacing(count=2)

            add_drag_float(name='specheat_1',default_value=1298,show=False,label='',format='%.0f',width=70)
            add_drag_float2(name='specheat_2',show=False,label='',format='%.0f',width=140)
            add_drag_float3(name='specheat_3',show=False,label='',format='%.0f',width=210)
            add_drag_float4(name='specheat_4',show=False,label='',format='%.0f',width=280)
            add_same_line(spacing=0)
            add_drag_float(name='specheat_ground',default_value=920,format='%.0f',show=True,label='Spec. Heat Cap. [J/kgºK]',width=70)
            
            add_spacing(count=2)

            add_drag_float(name='conduc_1',default_value=0.741,show=False,label='',format='%.3f',width=70)
            add_drag_float2(name='conduc_2',show=False,label='',format='%.3f',width=140)
            add_drag_float3(name='conduc_3',show=False,label='',format='%.3f',width=210)
            add_drag_float4(name='conduc_4',show=False,label='',format='%.3f',width=280)
            add_same_line(spacing=0)
            add_drag_float(name='conduc_ground',default_value=1.200,format='%.3f',show=True,label='Conductivity [W/mºK]',width=70)

            add_spacing(count=2)

            add_drag_float(name='thick_1',default_value=0.1524,show=False,label='',format='%.4f',width=70)
            add_drag_float2(name='thick_2',show=False,label='',format='%.4f',width=140)
            add_drag_float3(name='thick_3',show=False,label='',format='%.4f',width=210)
            add_drag_float4(name='thick_4',show=False,label='',format='%.4f',width=280)
            add_same_line(spacing=75)
            add_text('Layer Thickness [m]')

            add_spacing(count=2)

            add_drag_float(name='contactR_1',default_value=0.001,show=False,label='',format='%.3f',width=70)
            add_drag_float2(name='contactR_2',show=False,label='',format='%.3f',width=140)
            add_drag_float3(name='contactR_3',show=False,label='',format='%.3f',width=210)
            add_drag_float4(name='contactR_4',show=False,label='',format='%.3f',width=280)
            add_same_line(spacing=75)
            add_text('Therm Contact Resistance [0-1]')

            add_spacing(count=10)

            add_same_line(spacing=20)
            add_button('Previous##33',callback=hide3_3)
            add_same_line(spacing=440)
            add_button('Next##33',callback=show3_3,show=True)

        with tab('Thermal Profile Calculation Process', show=False):
            add_spacing(count=2)
            add_drag_int(name='numite',label='Max. number of recurring iterations',default_value=50,width=80,enabled=False)
            add_drag_int(name='convcrit',label='Convergence Criterion',default_value=50,width=80,enabled=False)

            add_spacing(count=5)

            add_text('Information about the recurring FTCS Finite Difference scheme:')
        
            add_spacing(count=5)

            add_text('Iteration Nº =')
            add_same_line(spacing=10)
            add_label_text(name='ite_number',label=' ')

            add_spacing(count=5)

            add_text('Frobenius Norm (convergence crit.) =')
            add_same_line(spacing=10)
            add_label_text(name='norm',label=' ')

            add_spacing(count=5)

            add_text('IMPORTANT EVENT :')
            add_same_line(spacing=10)
            add_label_text(name='event',label=' ')
            add_spacing(count=1)
            add_same_line(spacing=130)
            add_label_text(name='event2',label=' ')
            add_spacing(count=1)
            add_same_line(spacing=130)
            add_label_text(name='event3',label=' ')

            add_spacing(count=10)
            add_same_line(spacing=250)
            add_button('Recalculate',callback=show3_3,show=False)
            
            add_spacing(count=5)
            add_same_line(spacing=20)
            add_button('Previous##34',callback=hide3_4)
            add_same_line(spacing=440)
            add_button('Next##34',callback=show3_4,show=False)

        with tab('Thermal Profile / Time',show=False):
            add_plot('tempVStime',label='Temperature Evolution - NON-Ground Layers ONLY',height=360, \
                     x_axis_name='Time [h]',y_axis_name='Temperature [ºC]')

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##4',callback=hide4)
            add_same_line(spacing=186)
            add_button('Plot',callback=plot_thprofile)
            add_same_line(spacing=186)
            add_button('Next##4',callback=show4,show=False)

        with tab('Thermal Profile / Depth',show=False):
            add_spacing(count=2)

            add_combo('time_sel', \
                      label=' Time at which plotting the thermal profile [h]',width=100, \
                      callback=depthVStime)

            add_spacing(count=2)

            add_plot('tempVSdepth',label='Temperature Profile - ',height=320, \
                     x_axis_name='Temperature [ºC]',y_axis_name='Depth [m]')

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##5',callback=hide5)
            add_same_line(spacing=186)
            add_button('Export',callback=export,show=False)
            add_same_line(spacing=186)
            add_button('Next##5',callback=show5,show=False)

        with tab('Thermal Contraction Coeff',show=False):

            add_spacing(count=5)
        
            add_text('Linear coefficient of thermal contraction of the asphalt mix (1st layer only)')
        
            combo_item=['Linear coefficient input', \
                        'Mixture properties input', \
                        'Binder and Aggregates properties input']
        
            add_combo('thermal_coeff', items=combo_item, \
                      label='Input type',width=350, \
                      callback=coeff_sel)
        
            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)
        
            with collapsing_header('direct_in',label=combo_item[0],show=False,leaf=True):
        
                add_drag_float('bmix',default_value=2.351e-05, \
                                label='Linear coeff. of thermal contraction **Asphalt Mix** [1/ºC]', \
                                format='%.5e',width=150)
            
            with collapsing_header('partial_in',label=combo_item[1],show=False,leaf=True):
        
                add_drag_float('lbin##1',default_value=3.45e-4, \
                                label='Linear coeff. of thermal contraction **Binder** [1/ºK]', \
                                format='%.3e',width=150)
                add_drag_float('lag##1',default_value=9.5e-6, \
                                label='Linear coeff thermal contraction **Aggregates** [1/ºK]', \
                                format='%.3e',width=150)
                add_drag_float('vma',default_value=18.19, \
                                label='VMA intergranular void content in aggregate [%]',\
                                format='%.2f',width=150)
                add_drag_float('veff',default_value=81.71, \
                                label='Veff aggregate volume in mix [%]', \
                                format='%.2f',width=150)
            
            with collapsing_header('separate_in',label=combo_item[2],show=False,leaf=True):
        
                add_drag_float('lbin##2',default_value=3.45e-4, \
                                label='Linear coeff. of thermal contraction **Binder** [1/ºK]', \
                                format='%.3e',width=150)
                add_drag_float('lag##2',default_value=9.5e-6, \
                                label='Linear coeff thermal contraction **Aggregates** [1/ºK]', \
                                format='%.3e',width=150)
                add_drag_float('gsb',default_value=2.7, \
                                label='Aggregate bulk specific gravity', \
                                format='%.3e',width=150)
                add_drag_float('gb',default_value=1.03, \
                                label='Binder specific gravity', \
                                format='%.3e',width=150)
                add_drag_float('gmm',default_value=2.5, \
                                label='Asphalt mixture maximum specific gravity', \
                                format='%.3e',width=150)
                add_drag_float('av',default_value=7, \
                                label='Air volume [%]', \
                                format='%.1f',width=150)
                add_drag_float('bc',default_value=5, \
                                label='Binder content [%]', \
                                format='%.1f',width=150)
        
            add_spacing(count=4)
            add_drag_float(name='bmix_val',source='bmix_value',show=False, \
                           label='Coefficient used for calculation', \
                           format='%.5e',width=150,enabled=False) #no_input=True)
    
            add_spacing(count=10)
            add_same_line(spacing=20)
            add_button('Previous##51',callback=hide5_1)
            add_same_line(spacing=184)
            add_button('Validate',callback=valid_coeff)
            add_same_line(spacing=184)
            add_button('Next##51',callback=show5_1,show=False)

        with tab('Mechanical Properties Input',show=False):

            add_text('Select the material property you want to import:')

            mat_choice=['Creep Compliance','Relaxation Modulus']

            add_combo('mat_choice', items=mat_choice, \
                      label='Experimental data input',width=350,
                      callback=mat_input)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_table('mat_table',[],height=-50,show=False)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##6',callback=hide6)
            add_same_line(spacing=440)
            add_button('Next##6',callback=show6,show=False)

        with tab('Log-Log Optim Experimental Data',show=False):

            add_spacing(count=2)
            add_same_line(spacing=150)
            add_button('Plot Data', \
                       callback=plot_expe,callback_data=lambda:get_table_data('mat_table'))
            add_same_line(spacing=100)
            add_button('Plot Log-Log Data', \
                       callback=plot_expe,callback_data=lambda:get_table_data('mat_table'))
            add_spacing(count=2)

            add_plot('Experimental Data',height=325)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##7',callback=hide7)
            add_same_line(spacing=65)
            add_button('Optimization Poly(log-log)',callback=optim_expe, \
                       callback_data=lambda:get_table_data('mat_table'),show=False)
            add_same_line(spacing=65)
            add_button('Export##2',callback=optim_expe, \
                       callback_data=lambda:get_table_data('mat_table'),show=False)
            add_same_line(spacing=65)
            add_button('Next##7',callback=show7,show=False)

        with tab('Master Curve',show=False):

            add_spacing(count=2)
            add_same_line(spacing=80)
            add_button('Build Master Curve', enabled=False, \
                       callback=master_curve,callback_data=lambda:get_table_data('mat_table'))
            add_same_line(spacing=60)
            add_combo('tref##sel',width=80, \
                       label='Reference Temperature',callback=tref)
            add_spacing(count=2)

            add_plot('Master Curve##plot',xaxis_log_scale=True,yaxis_log_scale=True,height=325)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##8',callback=hide8)
            add_same_line(spacing=440)
            add_button('Next##8',callback=show8,show=False)
        
        with tab('Shift-Factor',show=False):

            add_spacing(count=2)
            add_same_line(spacing=80)
            # add_combo('shift_factor',width=260, items=['William-Landel-Ferry','Arrhenius'], \
            #            label='Shift-Factor Model',callback=shift,callback_data=lambda:get_table_data('mat_table'))
            add_combo('shift_factor',width=260, \
                       label='Shift-Factor Model',callback=shift,callback_data=lambda:get_table_data('mat_table'))
            add_spacing(count=2)

            add_plot('shift_plot',label='',height=275, \
                     x_axis_name='Temperature [ºC]',y_axis_name='log10(aT) [s]')

            add_spacing(count=2)
            add_same_line(spacing=100)
            add_drag_float2(name='wlf',source='coeff##wlf',show=False, \
                            label=' C1 and C2 [ºK] optimized values', \
                            format='%.2f',width=150,no_input=True)
            add_drag_float(name='arrh',source='coeff##arrh',show=False, \
                           label=' Material Activation Energy [J/mol]', \
                           format='%.2f',width=150,no_input=True)

            add_spacing(count=1)
            add_same_line(spacing=100)
            add_drag_float(name='r2',source='r_squared',show=False, \
                           label=' = R²', \
                           format='%.5f',width=150,no_input=True)

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##9',callback=hide9)
            add_same_line(spacing=440)
            add_button('Next##9',callback=show9,show=False)
        
        with tab('Prony Optimization##1',show=False):
            add_spacing(count=2)
            add_same_line(spacing=20)
            add_button('CCMC Optim', enabled=False, show=True, \
                       callback=ccmc_optim,callback_data=lambda:get_table_data('mat_table'))
            add_same_line(spacing=20)
            add_combo('nn##CCMC',width=50, items=['3','4','5','6'], show=True, \
                      label='General Mod. Power Law branch Nº',callback=ccmc_optim)
            add_same_line(spacing=15)
            add_button('Interconversion', show=False, \
                       callback=ccmc_optim,callback_data=lambda:get_table_data('mat_table'))
            add_spacing(count=2)

            add_plot('Material Optimization##plot',xaxis_log_scale=True,yaxis_log_scale=True,height=325, \
                     x_axis_name='Reduced Time [s]',y_axis_name='Creep Compliance [1/GPa]')

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##10',callback=hide10)
            add_same_line(spacing=440)
            add_button('Next##10',callback=show10,show=False)

        with tab('Prony Optimization##2',show=False):
            add_spacing(count=2)
            add_same_line(spacing=20)
            add_button('ErelMC Optim', enabled=False, show=True, \
                       callback=prony_optim2,callback_data=lambda:get_table_data('mat_table'))
            add_same_line(spacing=20)
            add_combo('nn##ErelMC',width=50, \
                      items=['4','5','6','7','8','9','10','11','12','13','14','15','16','18','20'], \
                      show=True,label='Generalized Maxwell Model branch number',callback=prony_optim2)
            add_spacing(count=2)

            add_plot('Material Optimization##plot2',xaxis_log_scale=True,yaxis_log_scale=True,height=325, \
                     x_axis_name='Reduced Time [s]',y_axis_name='Relaxation Modulus [GPa]')

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##11',callback=hide11)
            add_same_line(spacing=440)
            add_button('Next##11',callback=show11,show=False)

        with tab('Stress Profile / Time',show=False):
            add_plot('stressVStime',label='Stress Evolution - 1st layer ONLY',height=360, \
                     x_axis_name='Time [h]', \
                     y_axis_name='Stress [kPa]', \
                     #y_axis_name='Stress [PSI]', \
                     )

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##12',callback=hide12)
            add_same_line(spacing=186)
            add_button('Plot##2',callback=plot_stressprofile)
            add_same_line(spacing=186)
            add_button('Next##12',callback=show12,show=False)

        with tab('Stress Profile / Depth',show=False):
            add_spacing(count=2)

            add_combo('time_sel##2', \
                      label=' Time at which plotting the thermal profile [h]',width=100, \
                      callback=depthVSstress)

            add_spacing(count=2)

            add_plot('stressVSdepth',label='Stress Profile of 1st Layer - ',height=320, \
                     y_axis_name='Depth [m]', \
                     x_axis_name='Stress [kPa]', \
                     #x_axis_name='Stress [PSI]', \
                     )

            add_spacing(count=2)
            add_separator()
            add_spacing(count=2)

            add_spacing(count=1)
            add_same_line(spacing=20)
            add_button('Previous##13',callback=hide13)
            add_same_line(spacing=392)
            add_button('Export##3',callback=exportstress,show=False)

#show_documentation()
#show_logger()
#show_debug()
#show_metrics()
#show_about()

set_main_window_title('Asphalt Concrete Thermal Stress Calculator - ACTS Calc')
set_main_window_size(600,450)
set_theme('Light') # Classic; Light; Grey; Dark Grey; Dark; Dark 2

start_dearpygui(primary_window="ACTS")