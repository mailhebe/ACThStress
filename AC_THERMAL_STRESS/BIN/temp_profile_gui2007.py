import copy
import math
import cmath
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

####### dearpygui v0.6.415 #######
from dearpygui.core import *
from dearpygui.simple import *

####### for latest version #######
# from dearpygui.dearpygui import *

def wind_cubic(x,a,b,c,d):
    cubic_fit=a*(x**3)+b*(x**2)+c*x+d
    return cubic_fit

def wind_approx(days_to_compute,weather_data):

    wind_param=np.zeros((days_to_compute,4))
    wind_stderr=np.zeros((days_to_compute,4,4))

    for i in range(days_to_compute):
        wind_param[i,:],wind_stderr[i,:,:]=curve_fit(wind_cubic,weather_data[i]['Time [h]'], \
                                                    weather_data[i]['WindSpeed [mph]'])

    return wind_param,wind_stderr


def wind_poly_val(days_to_compute,weather_data,wind_param):
    wind_poly=np.zeros((days_to_compute,len(weather_data[0]['Time [h]'])))

    for i in range(days_to_compute):
        wind_poly[i,:] = wind_cubic(weather_data[i]['Time [h]'],*wind_param[i,:])

    for i in range(days_to_compute):
        for t in range(len(wind_poly[i,:])):
            if wind_poly[i,t]<0: wind_poly[i,t]=0

    return wind_poly


def r2_calc(days_to_compute,weather_data,wind_poly,wind_stderr):

    wind_r2=np.zeros((days_to_compute))

    for i in range(days_to_compute):
        residuals = weather_data[i]['WindSpeed [mph]']-wind_poly[i,:]
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((weather_data[i]['WindSpeed [mph]']-np.mean(weather_data[i]['WindSpeed [mph]']))**2)
        wind_r2[i] = 1 - (ss_res / ss_tot)

    return wind_r2

#############################################################################

def conversion_sewing(days_to_compute,weather_data,wind_poly):

    time = np.arange((len(weather_data[0]['Time [h]'])-1)*days_to_compute+1)

    temp_atm = [weather_data[0]['Tatm [ºC]'][0]]
    temp_dewpt = [weather_data[0]['DewPoint [ºF]'][0]]
    solar_rad = [weather_data[0]['SolarRad [W/m²]'][0]]
    wind_speed = [wind_poly.iloc[0,0]]

    for i in range(days_to_compute):
        temp_atm.extend(weather_data[i]['Tatm [ºC]'][1:].tolist())
        temp_dewpt.extend(weather_data[i]['DewPoint [ºF]'][1:].tolist())
        solar_rad.extend(weather_data[i]['SolarRad [W/m²]'][1:].tolist())
        wind_speed.extend(wind_poly.iloc[1:,i].tolist())
    
    temp_atm=np.array(temp_atm)
    temp_dewpt=np.array(temp_dewpt)
    solar_rad=np.array(solar_rad)
    wind_speed=np.array(wind_speed)

    temp_atm += 273.15 # [ºC] to [K]
    temp_dewpt = (temp_dewpt-32)/1.8 # [ºF] to [ºC]
    wind_speed *= (16/36) # [mph] to [m/s]

    return temp_atm,temp_dewpt,solar_rad,wind_speed

#############################################################################

def pavement_depth_const(pavement_thickness):
    
    pavement_depth=np.zeros(len(pavement_thickness))
    pavement_depth[0]=pavement_thickness[0]
    for k in range(1,len(pavement_thickness)):
        pavement_depth[k]=pavement_depth[k-1]+pavement_thickness[k]

    return pavement_depth

#############################################################################

def interior_nodes_cst(number_of_layers,pavement_depth,time_step,node_spacing,density,spec_heat_capacity,thermal_conductivity,interface_contact_res):

    delta_calc = np.zeros((number_of_layers))
    A = np.zeros((number_of_layers))
    B = np.zeros((number_of_layers))
    C=np.zeros((len(pavement_depth)))
    D=np.zeros((len(pavement_depth)))

    for layer in range(number_of_layers):
        delta_calc[layer] = (2*time_step)/(density[layer]*spec_heat_capacity[layer]*node_spacing)
        A[layer] = 1-(delta_calc[layer]*thermal_conductivity[layer]/node_spacing)
        B[layer] = 0.5*delta_calc[layer]*thermal_conductivity[layer]/node_spacing
    
    for interface in range(len(pavement_depth)):
        C[interface]=(2*node_spacing*thermal_conductivity[interface] \
                      + thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface]) \
                      / (2* (thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface] \
                           + node_spacing*thermal_conductivity[interface]+node_spacing*thermal_conductivity[interface+1]))
        D[interface]=(2*node_spacing*thermal_conductivity[interface+1] \
                      + thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface]) \
                      / (2* (thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface] \
                           + node_spacing*thermal_conductivity[interface]+node_spacing*thermal_conductivity[interface+1]))

    return delta_calc,A,B,C,D

def surface_nodes_cst(time_step,node_spacing,density,spec_heat_capacity,thermal_conductivity):

    delta_surf = (2*time_step)/(density[0]*spec_heat_capacity[0]*node_spacing)
    Bsurf = delta_surf*thermal_conductivity[0]/node_spacing
    return delta_surf,Bsurf

#############################################################################

def thermal_profile_calc(simul_prop,pavement_depth,ground_prop,surface_prop, \
    density,spec_heat_capacity,thermal_conductivity,interface_contact_res, \
    interp_import, \
    stefan_boltzmann,ref_temp_flow,ref_viscosity_flow,ref_conductivity_flow,ref_diffusivity_flow,ref_prandtl_flow, \
    numite, convcrit):

    number_of_layers = simul_prop[0]
    time_step = simul_prop[1]
    node_spacing = simul_prop[2]

    ground_layer_depth = ground_prop[0]
    deep_ground_temp = ground_prop[1]

    surface_albedo = surface_prop[0]
    surface_emissivity = surface_prop[1]
    sky_view_factor = surface_prop[2]
    solar_view_factor = surface_prop[3]
    characteristic_length = surface_prop[4]

    time_interp = interp_import.iloc[:,0].to_numpy()
    temp_atm_interp = interp_import.iloc[:,1].to_numpy()
    temp_dewpt_interp = interp_import.iloc[:,2].to_numpy()
    solar_rad_interp = interp_import.iloc[:,3].to_numpy()
    wind_speed_interp = interp_import.iloc[:,4].to_numpy()

    temp_sky_interp = temp_atm_interp*(0.004*temp_dewpt_interp+0.8)**(0.25)

    delta_surf = (2*time_step)/(density[0]*spec_heat_capacity[0]*node_spacing)
    Bsurf = delta_surf*thermal_conductivity[0]/node_spacing

    delta_calc = np.zeros((number_of_layers))
    A = np.zeros((number_of_layers))
    B = np.zeros((number_of_layers))
    C=np.zeros((len(pavement_depth)))
    D=np.zeros((len(pavement_depth)))

    for layer in range(number_of_layers):
        delta_calc[layer] = (2*time_step)/(density[layer]*spec_heat_capacity[layer]*node_spacing)
        A[layer] = 1-(delta_calc[layer]*thermal_conductivity[layer]/node_spacing)
        B[layer] = 0.5*delta_calc[layer]*thermal_conductivity[layer]/node_spacing
    
    for interface in range(len(pavement_depth)):
        C[interface]=(2*node_spacing*thermal_conductivity[interface] \
                      + thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface]) \
                      / (2* (thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface] \
                           + node_spacing*thermal_conductivity[interface]+node_spacing*thermal_conductivity[interface+1]))
        D[interface]=(2*node_spacing*thermal_conductivity[interface+1] \
                      + thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface]) \
                      / (2* (thermal_conductivity[interface]*thermal_conductivity[interface+1]*interface_contact_res[interface] \
                           + node_spacing*thermal_conductivity[interface]+node_spacing*thermal_conductivity[interface+1]))

    profile = np.zeros((len(time_interp),int(ground_layer_depth/node_spacing)+1))

    hrad = np.zeros((len(time_interp)))
    convective_heat_air = np.zeros((len(time_interp)))
    delta_t = np.zeros((len(time_interp)))
    Asurf = np.zeros((len(time_interp)))
    Csurf = np.zeros((len(time_interp)))

    normM = np.zeros((numite))

    unstab_flag = False
    convergence_flag = False

    for ite in range(numite):
        
        if ite==0:
            profile[-1,:].fill(deep_ground_temp)

        profile[:,-1].fill(deep_ground_temp)
        profile_old=copy.copy(profile)
        
        for time in range(0,profile.shape[0]):
            
            ##### hrad calculation #####
            hrad[time] = sky_view_factor*surface_emissivity*stefan_boltzmann*(profile[time-1,0]**2+temp_sky_interp[time]**2)*(profile[time-1,0]+temp_sky_interp[time])
            ############################
            
            ##### air film properties calculation #####
            film_temp = 0.5*(profile[time-1,0]+temp_atm_interp[time])            
            film_viscosity = ((film_temp-ref_temp_flow[0])/(ref_temp_flow[1]-ref_temp_flow[0])) \
                             *(ref_viscosity_flow[1]-ref_viscosity_flow[0])+ref_viscosity_flow[0]
            film_conductivity = ((film_temp-ref_temp_flow[0])/(ref_temp_flow[1]-ref_temp_flow[0])) \
                                *(ref_conductivity_flow[1]-ref_conductivity_flow[0])+ref_conductivity_flow[0]
            film_diffusivity = ((film_temp-ref_temp_flow[0])/(ref_temp_flow[1]-ref_temp_flow[0])) \
                                *(ref_diffusivity_flow[1]-ref_diffusivity_flow[0])+ref_diffusivity_flow[0]
            film_prandtl = ((film_temp-ref_temp_flow[0])/(ref_temp_flow[1]-ref_temp_flow[0])) \
                           *(ref_prandtl_flow[1]-ref_prandtl_flow[0])+ref_prandtl_flow[0]
            #film_beta = 1 / film_temp

            film_reynolds = wind_speed_interp[time] * characteristic_length / film_viscosity
            film_nusselt_laminar = 0.664 * film_reynolds**(0.5) * film_prandtl**(1/3)
            film_nusselt_turbulent = 0.037 * film_reynolds**(0.8) * film_prandtl**(1/3)

            if film_reynolds<5e5:
                film_nusselt_actual=film_nusselt_laminar
            else:
                film_nusselt_actual=film_nusselt_turbulent
                set_value('event','WATCH OUT : Turbulent flow!')
            ###########################################
            
            ##### hinf calculation (from air film prop) #####
            convective_heat_air[time] = film_nusselt_actual * film_conductivity / characteristic_length
            #################################################
            
            ##### CFL criterion validation #####
            delta_t[time]=(node_spacing**2)*density[0]*spec_heat_capacity[0] / (2*(hrad[time]*node_spacing+convective_heat_air[time]*node_spacing+thermal_conductivity[0]))
            stability = delta_t[time] > time_step
            
            if np.isnan(convective_heat_air[time])==True:
                print('NaN in h_inf calculation > originates from Nu calc... WHY?')
                print('Nu_actual =',film_nusselt_actual)
                print('film_temp =',film_temp)
                print('Tatm =',temp_atm_interp[time])

            if stability==False:
                unstab_flag = True
                break
            ####################################
            
            ##### A & C surface coefficients #####
            Asurf[time] = 1-(delta_surf*(hrad[time]+convective_heat_air[time]+thermal_conductivity[0]/node_spacing))
            Csurf[time] = delta_surf*(hrad[time]*temp_sky_interp[time]+convective_heat_air[time]*temp_atm_interp[time] + \
                               solar_view_factor*(1-surface_albedo)*solar_rad_interp[time])
            ######################################       
            
            for depth in range(profile.shape[1]):
                depth_m=round(depth*node_spacing,2)
                if depth_m==0:
                    profile[time,depth] = Asurf[time]*profile[time-1,depth]+Bsurf*profile[time-1,depth+1]+Csurf[time]
                    
                elif depth_m<=pavement_depth[-1]:
                    for interface in range(len(pavement_depth)):
                        if depth_m<pavement_depth[interface]:
                            profile[time,depth] = A[interface]*profile[time-1,depth] \
                                                + B[interface]*(profile[time-1,depth-1]+profile[time-1,depth+1])
                            break

                        elif depth_m==pavement_depth[interface]:
                            profile[time,depth] = C[interface]*profile[time-1,depth-1] \
                                                + D[interface]*profile[time-1,depth+1]
                            break

                elif depth_m!=ground_layer_depth:
                    profile[time,depth] = A[-1]*profile[time-1,depth] \
                                        + B[-1]*(profile[time-1,depth-1]+profile[time-1,depth+1])

                else:
                    profile[time,depth] = deep_ground_temp

        if unstab_flag == True:
            set_value('event','UNSTABLE FD scheme, breaking!')
            set_value('event2','Change TIME and/or SPACE discretization parameter(s).')
            set_value('event3','See CFL [Courant-Friedrichs-Lewy] criterion for info.')
            break
        
        profile_diff=profile_old-profile
        normM[ite] = np.linalg.norm(profile_diff,'fro')

        set_value('ite_number',str(ite+1))
        set_value('norm',str(normM[ite]))
        
        if normM[ite]<=convcrit:
            convergence_flag = True
            set_value('event','CONVERGENCE CRITERION Fnorm(X)<50 REACHED!')
            break

    return profile,normM,ite,unstab_flag,convergence_flag