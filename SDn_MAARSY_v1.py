# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 23:31:57 2023

@author: ma042
"""

import matplotlib.pyplot as plt
import h5py as h5
import numpy as np
import pandas as pd

import os
import datetime
from datetime import datetime, timedelta
from matplotlib.dates import date2num
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
import imageio


## self-generated codes
# from extract_winds import extract_winds
from structure_functions import calculate_structure_function_q, calculate_structure_function_q_loc, calculate_structure_function_q_loc_binned

################################################################################
### FUNCTIONS
################################################################################
def second_order_fit(r, a, b, c):
        return a*r**(2/3) + b*r**2 - c*r**2*np.log(r)
        # return b*r**2 - c*r**2*np.log(r)

def cube(x):
    if x >= 0:
        return x**(1/3)
    elif x < 0:
        return -(abs(x)**(1/3))




# Read the dates from the file
time_sequence_15min = []
with open('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/date_list_15min_2019-2022.txt', 'r') as file:
    for line in file:
        # Parse each line as a datetime object
        time = datetime.strptime(line.strip(), '%Y-%m-%d %H:%M:%S')
        time_sequence_15min.append(time)

# for timestamp in time_sequence_15min[:100]:  # Display the first 100 timestamps
#     print(timestamp)



W_15min = np.load('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/U0mod_15min_2019-2022.npy')

H = np.load('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/height.npy')








# #####################################################################################################
# # Plot
# plt.figure(figsize=(10, 6))
# plt.pcolormesh(WINDS.T, shading='auto')

# plt.colorbar(label='Wind Velocity Modulus (m/s)')
# plt.title('Wind Velocity Modulus')
# plt.gca().xaxis_date()  # Set the x-axis as dates
# # plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
# plt.tight_layout()
# plt.show()


 
# ## TURBULENCE INTENSITY OF THE WINDS for every altitude
# TI = np.nanstd(WINDS,axis=0)/np.nanmean(WINDS,axis=0)*100

# # plt.figure()
# # plt.plot(H,TI,'.')
# # plt.xlim(0,20)
# # plt.ylim(0, 100)
# # plt.grid()


#####################################################################################################
# SET PARAMETERS:
q = 2  # Replace with the desired q value
# max_tau = 1000  # Replace with the maximum value of tau you want to calculate
# max_tau = 400  # Replace with the maximum value of tau you want to calculate
# res = 15
max_tau = 200  # Replace with the maximum value of tau you want to calculate
res = 15
tau = np.arange(res*60,(max_tau+1)*res*60,res*60) # in seconds
i_tt = np.where((tau == 1 * 60 * 60) | (tau == 24 * 60 * 60) | (tau == 3 * 24 * 60 * 60) | (tau == 30 * 24 * 60 * 60))

# values for binned taylor approximation
r_min = 10 * 1000 # minimum value in meters. It is constrained by MAARSY measurement area.
r_max = 10000 * 1000
num_bins = 500
bin_edges = np.logspace(np.log10(r_min), np.log10(r_max), num_bins + 1) 


# #####################################################################################################

# #####################################################################################################
# ### OLD: Functions that call to SFs estimation functions and plot results (not currently run)
# #####################################################################################################

def temp_SFs():
    # TEMPORAL SFs
    ## Altitude range: These ranges allow to clearly distinguish different behaviors
    ## in termporal third-order SFs.
    ih0 = np.squeeze(np.where((H>=6) & (H<6.5)))
    ih1 = np.squeeze(np.where((H>=6.5) & (H<11)))
    ih2 = np.squeeze(np.where((H>=6) & (H<10)))
    ihT = np.squeeze(np.where((H>=1) & (H<20)))
    
    for ih in [ih2]:
        print(ih)
        fig, ax = plt.subplots()
        Dq_13 = []
        DQ_loc = []
        r = []
        for c in ih:
            print(c)
            wind_speed_data = W_15min[:, c]
            D_q = calculate_structure_function_q(wind_speed_data, q, max_tau)
           
            d_q13 = []
            for d_q in D_q:
                d_q13.append(cube(d_q))
            Dq_13.append(np.array(d_q13))
    
            DQ_loc.append(Dq_13)
    
            # plt.semilogx(tau, d_q13, '.', label='', alpha=0.1)
        
        # Calculate mean and standard deviation
        mean_Dq_13 = np.nanmean(np.vstack(Dq_13), axis=0)
        median_Dq_13 = np.nanmedian(np.vstack(Dq_13), axis=0)
        std_Dq_13 = np.nanstd(np.vstack(Dq_13), axis=0)
        
        # Plot mean and fill_between
        plt.semilogx(tau, mean_Dq_13, '.k', label='Mean')
        plt.semilogx(tau, median_Dq_13, 'xr', label='Median')
        plt.fill_between(tau, mean_Dq_13 + std_Dq_13, mean_Dq_13 - std_Dq_13, color='gray', alpha=0.3, label='± 1 std')
        
        plt.ylim([-10, 10])
        plt.legend()
        plt.grid()
        
        plt.xlabel(r'$\tau$ (s)')
        plt.ylabel(r'$D_3^{1/3}$ (m/s)')
        plt.title('Heights: %1.1f-%1.1f km'%(min(H[ih]),max(H[ih])))
    
    
        # Add vertical lines and text annotations
        for idx in i_tt[0]:
            plt.axvline(x=tau[idx], color='k', linestyle='--', linewidth=2.0)
            plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d', '30 d'][i_tt[0].tolist().index(idx)], color='k', rotation=90,
                      verticalalignment='bottom', horizontalalignment='right')
    
        plt.show()
    # plt.savefig("C:/Poblet/IAP/work/N-order_sf/figuritas/D3_tau_MAARSY-10_13km_%4d_maxtau%d.jpg"%(y,max_tau), dpi=400)


def third_order_spatial_SFs():
    # Third-order structure functions
    ### Instead of the averaged result, particular heights are used
    for ih in range(len(ihs)):
        r2Ds = r2D[ih]
        r2Dcounts = r2Dcount[ih]
        DQ2Ds = DQ2D[ih]
        
        # Positive - Negative discrimination - 1D curves
        Dq_N = np.empty((len(bin_edges)))*np.nan
        Dq_P = np.empty((len(bin_edges)))*np.nan
        
        r = np.nanmean(r2Ds, axis=1)
        DQ_loc = np.nanmean(DQ2Ds, axis=1)
        
        iN = np.where(DQ_loc < 0)
        iP = np.where(DQ_loc >=0)
        
        Dq_N[iN] = DQ_loc[iN]
        Dq_P[iP] = DQ_loc[iP]
        
        fig, ax = plt.subplots()
        
        # plt.loglog(r/1000, (r/1000)**(2/3), '--k', linewidth=2, label='$\sim s^{2/3}$')
        # plt.loglog(r/1000, 0.04*(r/1000)**3, '--k', linewidth=2, label='$\sim s^3$')
        
        plt.loglog(r/1000, np.abs(Dq_N), '.r', linewidth=2, label = r'$D_3<0$')
        plt.loglog(r/1000, np.abs(Dq_P), '.b', linewidth=1, label = r'$D_3\geq 0$')
        
        # plt.loglog(r/1000, np.abs(DQ2D), '.', color = 'gray',alpha = 0.5)
     
        ## limit the s range to do the fit
        is_fit = np.where(bin_edges <= 5000*1000)
        # popt, pcov = curve_fit(second_order_fit, bin_edges[is_fit]/1000, np.abs(Dq_P[is_fit])) # exclude last point because they weight too much on the fit
        # print(popt)
        # plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**(2/3) + popt[1]*(bin_edges/1000)**2 - popt[2]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), color='orange', linewidth=2, label=r'$\hat{D}_2 = a_0 s^{2/3} + a_1 s^2 - a_2 s^2 \ln{(s)}$')    
      
     
        plt.xlabel('$s$ (km)')
        plt.ylabel('(m$^3$s$^{-3}$)')
        # plt.title('$D_3 =<(u_r´-u_r)^3>$. H: %1.1f-%1.1f km'%(H[ihs[0]],H[ihs[-1]]))
        plt.title('$D_3 =<(u_r´-u_r)^2>$. H: %1.1f km'%(H[ihs[ih]]))
        
        
        plt.xlim(1,20000)
        plt.ylim(1e-2, 1e3)
        
        plt.grid()
        plt.show()
        plt.legend()
        # plt.savefig("C:/Poblet/IAP/work/N-order_sf/figuritas/v_1/single_heights_2019-2022/D2/D3_s_MAARSY%02d.jpg"%(ih), dpi=400)
    # plt.close('all')        
            

#####################################################################################################
# Call functions that calculate SFs
#####################################################################################################

# ihs =  np.squeeze(np.where((H>=5) & (H<15)))
ihs =  np.squeeze(np.where((H>=8.5) & (H<12.1)))
q = 2 # order of the SFs

row = 0
for ih in [ihs]:    
    print(ih)
    r2D_list = []
    r2Dcount_list = []
    Vmn2Dcount_list = []
    DQ2D_loc_list = []
    for c in ih:
        print(c)
        wind_speed_data = W_15min[:, c]
        Dq_loc_binned, r_loc_binned, rcount_loc_binned, vmncount_loc_binned = calculate_structure_function_q_loc_binned(wind_speed_data, 
                                                                                                   q, 
                                                                                                   max_tau, 
                                                                                                   tau, 
                                                                                                   r_min, 
                                                                                                   r_max, 
                                                                                                   num_bins, 
                                                                                                   central_tendency='mean')
        
        r2Dcount_list.append(rcount_loc_binned)   
        Vmn2Dcount_list.append(np.array(vmncount_loc_binned))   

        DQ2D_loc_list.append(Dq_loc_binned)
        r2D_list.append(r_loc_binned)
 
## Calculate mean 2D arrays using different altitudes

### Convert the list of arrays into a single numpy array with one additional dimension
# Vmn2Dcount = np.stack(Vmn2Dcount_list)
r2Dcount = np.stack(r2Dcount_list)

r2D = np.stack(r2D_list)    
DQ2D = np.stack(DQ2D_loc_list)

# DQ2D = np.nanmean(DQ2D,axis=0)

# #############################################################################
# # Analysis of advection velocitites at particular lags

# ## Complete signal for some taus
# plt.figure()
# plt.plot(Vmn2Dcount_list[1][5],'k.')
# plt.plot(Vmn2Dcount_list[1][200],'r.')
# plt.plot(Vmn2Dcount_list[1][350],'b.')
# plt.grid()

# ## mean and median for all taus
# mean_Vnm = np.empty(len(Vmn2Dcount_list[0]))*np.nan
# median_Vnm = np.empty(len(Vmn2Dcount_list[0]))*np.nan

# for im in range(len(Vmn2Dcount_list[0])):
#     mean_Vnm[im] = np.nanmean(Vmn2Dcount_list[0][im])
#     median_Vnm[im] = np.nanmedian(Vmn2Dcount_list[0][im])
    

    
# plt.figure()
# plt.plot(tau/60/60, mean_Vnm,'.r', label = 'mean')
# plt.plot(tau/60/60, median_Vnm,'.b', label = 'median')
# plt.axhline(y=np.nanmedian(wind_speed_data), color='k', linestyle='--', label = 'median (global)')
# plt.axhline(y=np.nanmean(wind_speed_data), color='k', linestyle='-', label = 'mean (global)')

# plt.xlabel(r'$\tau$ (h)' )

# plt.legend()

# ## distributions for particular taus
# plt.figure()
# L = [r'$\tau=$1 h', r'$\tau=$1 d', r'$\tau=$3 d']
# c = 0
# for idx in i_tt[0]:

#     pdf_norm, bins = np.histogram(Vmn2Dcount_list[1][idx], bins=np.linspace(0, 70+1, 500), density=False)

#     std_norm_plot = np.nanstd(Vmn2Dcount_list[1][idx])

#     # mean_Da_norm_plot = np.nanmean(Da_norm,axis=0)/std_Da_norm_plot
#     # median_Da_norm_plot = np.nanmedian(Da_norm)/std_Da_norm_plot

#     plt.step(bins[:-1]/std_norm_plot, pdf_norm, linewidth=2, label=L[c])

#     # plt.step(bin_edges/1000, r2Dcount[10,:,idx], label=L[c])
#     # plt.plot(bins_s[1:], pdf_s)
    
#     c += 1
# plt.title(r'histogram of $V_{m,n}$ for fixed $\tau$. H: %1.1f km'%(H[ihs[1]]))
# plt.xlabel(r'$V_{m,n}/\sigma_{V_{m,n}}$')
# # plt.xscale('log')
# plt.legend()
# plt.grid()


# # tau vs s
# plt.figure()
# plt.pcolormesh(tau, bin_edges/1000, r2Dcount[1], shading='auto', cmap='plasma')
# plt.colorbar()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\tau$ (s)')
# plt.ylabel(r'$s$ (km)')
# plt.title(r'# of values for each $(\tau,s)$ bin')







#############################################################################
# Second-order SFs plots

# A1: Lindborg (2018) - COMPARISON WITH AVERAGED CURVE
df = pd.read_csv('C:/Poblet/IAP/work/N-order_sf/DllplutDtt_lindborg_2018.csv')
dffs_l = pd.read_csv('C:/Poblet/IAP/work/N-order_sf/Dll_AMDAR_fig3_FREHLICH2010.csv')
dffs_t = pd.read_csv('C:/Poblet/IAP/work/N-order_sf/Dtt_AMDAR_fig3_FREHLICH2010.csv')

# Use the .to_numpy() method to convert the selected columns to a numpy array
# s = df[1].to_numpy()
D2strat = df.to_numpy()
D2fs_l = dffs_l.to_numpy()
D2fs_t = dffs_t.to_numpy()

# Create an empty list to store filenames of saved images
filenames = []

# Your loop with modifications to save figures
for ih in range(len(ihs)):
    r2Ds = r2D[ih]
    r2Dcounts = r2Dcount[ih]
    DQ2Ds = DQ2D[ih]
    
    ## Calculate temporal SFs to compare in the plot
    wind_speed_data = W_15min[:, ih]
    D_qt = calculate_structure_function_q(wind_speed_data, q, max_tau) 
   
    # What central tendency value of the r distribution performs better?
    r = np.nanmean(r2Ds, axis =1) # mean
    # r = np.nanmedian(r2Ds, axis =1) # median
    
    DQ_loc = np.nanmean(DQ2Ds, axis =1)    

    fig, ax = plt.subplots()
   
    is_fit = np.where(bin_edges <= 5000*1000)
    popt, pcov = curve_fit(second_order_fit, bin_edges[is_fit]/1000, np.abs(DQ_loc[is_fit]))

    # full SFs
    plt.loglog(tau*np.nanmedian(wind_speed_data)/1000, D_qt, '.k',  label = '$D_2$')
    plt.loglog(r/1000, DQ_loc, '.b', label = '$D^{(l)}_2$')
    # plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**(2/3) + popt[1]*(bin_edges/1000)**2 - popt[2]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), color='orange', linewidth=2, label=r'$\hat{D}_2 = a_0 s^{2/3} + a_1 s^2 - a_2 s^2 \ln{(s)}$')    
 
    # # mesoscale term
    # plt.loglog(r/1000, DQ_loc - popt[1]*(bin_edges/1000)**2 - popt[2]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), '.b', linewidth=1, label = '$D_2 - a_1 s^2 - a_2 s^2 \ln{(s)}$')
    # plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**(2/3), color='orange', linewidth=2, label=r'$\hat{D}_2 = a_0 s^{2/3}$')    

    # # Dll+Dtt in the stratosphere    
    # plt.loglog(D2strat[:,0],D2strat[:,1],'k') # 
    plt.loglog(D2strat[:,0],D2strat[:,1]/2, linewidth=3, label=r'$(D_{ll}+D_{tt})/2$ from Lindborg (2014)') # Divide by 2 to get an average of the two components
    plt.loglog(D2fs_l[:,0],D2fs_l[:,1], linewidth=3, label=r'$D_{ll}$ from Frehlich & Sharman (2010)') # 
    plt.loglog(D2fs_t[:,0],D2fs_t[:,1], linewidth=3, label=r'$D_{tt}$ from Frehlich & Sharman (2010)') # 
    
    
    plt.xlabel('$s$ (km)')
    plt.ylabel('(m$^2$s$^{-2}$)')
    plt.title('$D_2 =<(u_r´-u_r)^2>$. H: %1.1f km'%(H[ihs[ih]]))
    plt.xlim(1,20000)
    # plt.ylim(1e-0, 1e2)
    plt.ylim(1e-1, 1e4)

    plt.grid()
    plt.legend()
    
#     # Save the current figure
#     # filename = "C:/Poblet/IAP/work/N-order_sf/figuritas/v_1/single_heights_2019-2022/D2/mean/D2_s_MAARSY{:02d}.jpg".format(ih)
#     filename = "C:/Poblet/IAP/work/N-order_sf/figuritas/v_2/v_2_15min_res_2019-2022/without_fits/D2_s_MAARSY{:02d}.jpg".format(ih)
#     plt.savefig(filename, dpi=400)
#     filenames.append(filename)

#     # Close the figure to avoid displaying all figures
#     plt.close()


# # Create a GIF from saved images
# with imageio.get_writer('C:/Poblet/IAP/work/N-order_sf/figuritas/v_2/v_2_15min_res_2019-2022/without_fits/figures.gif', mode='I', duration=0.5) as writer:
#     for filename in filenames:
#         image = imageio.imread(filename)
#         writer.append_data(image)

























# ## PLOTS
# # ## D3
# plt.figure(figsize=(10, 6))
# plt.pcolormesh(tau, bin_edges/1000, DQ2D, shading='auto', cmap='seismic_r', vmin=-1000, vmax=1000)
# plt.colorbar(label=r'$D_3(\tau,s)$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\tau$ (s)')
# plt.ylabel(r'$s$ (km)')
# plt.title(r'Third-order SF $D_3(\tau,s)$')
# # Add vertical lines and text annotations
# for idx in i_tt[0]:
#     plt.axvline(x=tau[idx], color='orange', linestyle='--', linewidth=2.0)
#     plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d', '30 d'][i_tt[0].tolist().index(idx)], color='k', rotation=90,
#               verticalalignment='top', horizontalalignment='right')

# # plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
# plt.show()

# # ## R vs tau : # of values
# plt.figure(figsize=(10, 6))
# plt.pcolormesh(tau, bin_edges/1000, r2Dcount, shading='auto', cmap='plasma')
# plt.colorbar()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\tau$ (s)')
# plt.ylabel(r'$s$ (km)')
# plt.title(r'# of values for each $(\tau,s)$ bin')

# # plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
# # plt.tight_layout()

# # Add vertical lines and text annotations
# for idx in i_tt[0]:
#     plt.axvline(x=tau[idx], color='orange', linestyle='--', linewidth=2.0)
#     plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d'][i_tt[0].tolist().index(idx)], color='orange', rotation=90,
#               verticalalignment='bottom', horizontalalignment='right')
# plt.show()


# # ## single tau plots
# plt.figure(figsize=(10, 6))
# L = [r'$\tau=$1 h', r'$\tau=$1 d', r'$\tau=$3 d']
# c = 0 
# for idx in i_tt[0]:
#     # ## Central tendendy estimates of s distributions for fixed taus    
#     # pdf_s, bins_s = np.histogram(r2Dcount[:,idx], bins=bin_edges/1000, density=False)
  
    
    
#     plt.step(bin_edges/1000, r2Dcount[10,:,idx], label=L[c])
#     # plt.plot(bins_s[1:], pdf_s)
    
#     c += 1
# plt.title(r'histogram of $s$ for fixed $\tau$')
# plt.xlabel(r'$s$ (km)')
# plt.xscale('log')
# plt.legend()
# plt.grid()





# # ### Kurtosis
# # plt.figure()
# # for c in [15,16,17,18]:
# #     wind_speed_data = WINDS[:,c]
    
# #     # Example usage:
# #     q = 3  # Replace with the desired q value
# #     max_tau = 288*15  # Replace with the maximum value of tau you want to calculate
# #     D_2 = calculate_structure_function_q(wind_speed_data, 2, max_tau)
# #     D_4 = calculate_structure_function_q(wind_speed_data, 4, max_tau)

    
    
    
# #     # D_q_loc = calculate_structure_function_q_loc(wind_speed_data, q, max_tau)
    
# #     # D_q contains D_q(1) to D_q(max_tau)
# #     # D_q_loc contains D_q^loc(r) for each tau from 1 to max_tau
    
    
# #     # plt.figure()
# #     # for d in D_q_loc:
# #     #     plt.plot(d,'-o')
    
# #     tau = np.arange(5*60,(max_tau+1)*5*60,5*60) 
    
    
# #     plt.semilogx(tau,D_4/D_2**2,'.')
    
# # plt.grid()


