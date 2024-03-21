# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 23:31:57 2023

@author: ma042
"""

import matplotlib.pyplot as plt
import h5py as h5
import numpy as np

import os
import datetime
from datetime import datetime
from matplotlib.dates import date2num
from scipy.optimize import curve_fit

## self-generated codes
from extract_winds import extract_winds
from structure_functions import calculate_structure_function_q, calculate_structure_function_q_loc, calculate_structure_function_q_loc_binned

################################################################################
### FUNCTIONS
################################################################################
def second_order_fit(r, b,c):
        # return a*r**(2/3) + b*r**2 - c*r**2*np.log(r)
        return b*r**2 - c*r**2*np.log(r)

def cube(x):
    if x >= 0:
        return x**(1/3)
    elif x < 0:
        return -(abs(x)**(1/3))


year = [2019, 2020, 2021, 2022]
months = [1,2,3,4,5,6,7,8,9,10,11,12]
# dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]

# for y in [2019, 2020, 2021,2022]:
#     year = [y]
#     months = [1,2,3,4,5,6,7,8,9,10,11,12]
#     dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]
    
# # SUMMERS
# year = [2019, 2020, 2021, 2022]
# months = [6,7,8]

# # WINTERS
# year = [2019, 2020, 2021, 2022]
# months = [11,12,1]



path = 'C:/Poblet/IAP/work/N-order_sf/MAARSY/'

U0, V0, WINDS, H, date_list =  extract_winds(ipath = path,
                                              year = year,
                                              month=months,                                              
                                              component=2) # 0: zonal, 1:meridional, 2:horizontal absolute value, 3: full amplitude


#####################################################################################################
# # SAVE FILES
# with open('C:/Poblet/IAP/work/N-order_sf/data/date_list_2022.txt', 'w') as f:
#     f.write('\n'.join([date.strftime('%Y-%m-%d %H:%M:%S') for date in date_list]))


# np.save('C:/Poblet/IAP/work/N-order_sf/data/WIND2022.npy', WINDS)
# np.save('C:/Poblet/IAP/work/N-order_sf/data/height.npy', H)


# # LOAD FILES
# date_list = [line.strip() for line in open('C:/Poblet/IAP/work/N-order_sf/data/date_list_2022.txt')]
# WINDS = np.load('C:/Poblet/IAP/work/N-order_sf/data/WIND2022.npy')
# H = np.load('C:/Poblet/IAP/work/N-order_sf/data/height.npy')

#####################################################################################################



## Checking winds
date_nums = date2num(date_list)  # Convert to matplotlib date format

# Plot
# plt.figure(figsize=(10, 6))
# plt.pcolormesh(date_nums, H, WINDS.T, shading='auto')
# plt.colorbar(label='Wind Velocity Modulus (m/s)')
# plt.xlabel('Date')
# plt.ylabel('Height (km)')
# plt.title('Wind Velocity Modulus')
# plt.gca().xaxis_date()  # Set the x-axis as dates
# # plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
# plt.tight_layout()
# plt.show()


 
## TURBULENCE INTENSITY OF THE WINDS for every altitude

TI = np.nanstd(WINDS,axis=0)/np.nanmean(WINDS,axis=0)*100

# plt.figure()
# plt.plot(H,TI,'.')
# plt.xlim(0,20)
# plt.ylim(0, 100)
# plt.grid()


#####################################################################################################
# SET PARAMETERS:
q = 3  # Replace with the desired q value
# max_tau = 1000  # Replace with the maximum value of tau you want to calculate
max_tau = 400  # Replace with the maximum value of tau you want to calculate
tau = np.arange(15*60,(max_tau+1)*15*60,15*60) # in seconds
i_tt = np.where((tau == 1 * 60 * 60) | (tau == 24 * 60 * 60) | (tau == 3 * 24 * 60 * 60) | (tau == 30 * 24 * 60 * 60))

# values for binned taylor approximation
r_min = 10 * 1000 # minimum value in meters. It is constrained by MAARSY measurement area.
r_max = 10000 * 1000
num_bins = 500
bin_edges = np.logspace(np.log10(r_min), np.log10(r_max), num_bins + 1) 

# #####################################################################################################

# #####################################################################################################
# ### Third order SFs plots
# #####################################################################################################
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
        wind_speed_data = WINDS[:, c]
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

#####################################################################################################

#####################################################################################################
# SPATIAL THIRD ORDER SFs
ihs =  np.squeeze(np.where((H>=5) & (H<15)))
# ihs =  np.squeeze(np.where((H>=1) & (H<6.5)))
q = 3

row = 0
for ih in [ihs]:    
    print(ih)
    r2D_list = []
    r2Dcount_list = []
    DQ2D_loc_list = []
    for c in ih:
        print(c)
        wind_speed_data = WINDS[:, c]
        Dq_loc_binned, r_loc_binned, rcount_loc_binned = calculate_structure_function_q_loc_binned(wind_speed_data, 
                                                                                                   q, 
                                                                                                   max_tau, 
                                                                                                   tau, 
                                                                                                   r_min, 
                                                                                                   r_max, 
                                                                                                   num_bins, 
                                                                                                   median=False)
        
        DQ2D_loc_list.append(Dq_loc_binned)
        r2D_list.append(r_loc_binned)
        r2Dcount_list.append(rcount_loc_binned)   
 
## Calculate mean 2D arrays using different altitudes

### Convert the list of arrays into a single numpy array with one additional dimension
r2D = np.stack(r2D_list)    
r2Dcount = np.stack(r2Dcount_list)
DQ2D = np.stack(DQ2D_loc_list)

### Calculate the mean along the added dimension (axis=0)
r2D = np.nanmean(r2D,axis=0)
r2Dcount = np.nansum(r2Dcount,axis=0)
DQ2D = np.nanmean(DQ2D,axis=0)

### Instead of the averaged result, particular heights are used

for ih in range(len(ihs)):
    r2Ds = r2D[ih]
    r2Dcounts = r2Dcount[ih]
    DQ2Ds = DQ2D[ih]
    
    # Positive - Negative discrimination - 1D curves
    Dq_N = np.empty((len(bin_edges)))*np.nan
    Dq_P = np.empty((len(bin_edges)))*np.nan
    
    r = np.nanmean(r2Ds, axis =1)
    DQ_loc = np.nanmean(DQ2Ds, axis =1)
    
    iN = np.where(DQ_loc < 0)
    iP = np.where(DQ_loc >=0)
    
    Dq_N[iN] = DQ_loc[iN]
    Dq_P[iP] = DQ_loc[iP]
    
    
    fig, ax = plt.subplots()
    plt.loglog(r/1000, 0.03*(r/1000), '--k', linewidth=2, label='$\sim s^3$')
    # plt.loglog(r/1000, 0.04*(r/1000)**3, '--k', linewidth=2, label='$\sim s^3$')
    
    plt.loglog(r/1000, np.abs(Dq_N), '.r', linewidth=2, label = '$D_3<0$')
    plt.loglog(r/1000, np.abs(Dq_P), '.b', linewidth=2, label = '$D_3>=0$')
    
    # plt.loglog(r/1000, np.abs(DQ2D), '.', color = 'gray',alpha = 0.5)
    
    plt.xlabel('$s$ (km)')
    plt.ylabel('(m$^3$s$^{-3}$)')
    # plt.title('$D_3 =<(u_r´-u_r)^3>$. H: %1.1f-%1.1f km'%(H[ihs[0]],H[ihs[-1]]))
    plt.title('$D_3 =<(u_r´-u_r)^3>$. H: %1.1f km'%(H[ihs[ih]]))
    
    
    plt.xlim(1,20000)
    plt.ylim(1e-2, 1e3)
    
    plt.grid()
    plt.show()
    # plt.legend()
    plt.savefig("C:/Poblet/IAP/work/N-order_sf/figuritas/v_1/single_heights_2019/D3_s_MAARSY%02d.jpg"%(ih), dpi=400)
        
            
       
## PLOTS
# ## D3
plt.figure(figsize=(10, 6))
plt.pcolormesh(tau, bin_edges/1000, DQ2D, shading='auto', cmap='seismic_r', vmin=-1000, vmax=1000)
plt.colorbar(label=r'$D_3(\tau,s)$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\tau$ (s)')
plt.ylabel(r'$s$ (km)')
plt.title(r'Third-order SF $D_3(\tau,s)$')
# Add vertical lines and text annotations
for idx in i_tt[0]:
    plt.axvline(x=tau[idx], color='orange', linestyle='--', linewidth=2.0)
    plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d', '30 d'][i_tt[0].tolist().index(idx)], color='k', rotation=90,
              verticalalignment='top', horizontalalignment='right')

# plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
plt.show()

# ## R vs tau : # of values
plt.figure(figsize=(10, 6))
plt.pcolormesh(tau, bin_edges/1000, r2Dcount, shading='auto', cmap='plasma')
plt.colorbar()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\tau$ (s)')
plt.ylabel(r'$s$ (km)')
plt.title(r'# of values for each $(\tau,s)$ bin')

# plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
# plt.tight_layout()

# Add vertical lines and text annotations
for idx in i_tt[0]:
    plt.axvline(x=tau[idx], color='orange', linestyle='--', linewidth=2.0)
    plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d'][i_tt[0].tolist().index(idx)], color='orange', rotation=90,
              verticalalignment='bottom', horizontalalignment='right')
plt.show()


# ## single tau plots
plt.figure(figsize=(10, 6))
L = [r'$\tau=$1 h', r'$\tau=$1 d', r'$\tau=$3 d']
c = 0 
for idx in i_tt[0]:
    plt.plot(bin_edges/1000, r2Dcount[:,idx], label=L[c])
    # plt.axvline(x=np.nanmax(r2Dcount[:,idx]), color='k', linestyle='--', linewidth=2.0)
    c += 1
plt.xlabel(r'$s$ (km)')
plt.xscale('log')
plt.legend()
plt.grid()


# ## Determine the slope of R vs tau plot
mean_r2Dcount = []
median_r2Dcount = []
max_r2Dcount = []
for idy in range(0,len(tau)):
    mean_r2Dcount.append(np.nanmean(r2Dcount[:,idy]))
    median_r2Dcount.append(np.nanmedian(r2Dcount[:,idy]))
    max_r2Dcount.append(np.nanmax(r2Dcount[:,idy]))
plt.figure(figsize=(10, 6))
# plt.plot(tau,np.array(mean_r2Dcount),'o')
# plt.plot(tau,np.array(median_r2Dcount),'o')
plt.plot(tau,np.array(max_r2Dcount),'o')
# plt.xscale('log')
# plt.yscale('log')


  

    


















    


# ################################################################################
# # SPATIAL SECOND ORDER SFs
# ihs =  np.squeeze(np.where((H>=11) & (H<13)))
# # ihs =  np.squeeze(np.where((H>=1) & (H<6.5)))
# q = 2

# DQ = np.empty((len(ihs),len(bin_edges)))*np.nan
# row = 0
# for ih in [ihs]:    
#     print(ih)

#     dq_loc_binned = np.empty(len(bin_edges))*np.nan
#     DQ_loc = []

#     r = []
#     for c in ih:
#         print(c)
#         wind_speed_data = WINDS[:, c]
#         Dq_loc_mean, r_loc_mean = calculate_structure_function_q_loc(wind_speed_data, 
#                                                                                                   q, 
#                                                                                                   max_tau,
#                                                                                                   tau)
        
#         # With which array I am working
#         dtw = Dq_loc_mean
#         r_l = r_loc_mean

#         DQ_loc.append(dtw)
#         r.append(r_l)




# fig, ax = plt.subplots()
# for p in range(len(r)):
    
#     # Positive - Negative discrimination
#     Dq_N = np.empty((len(DQ_loc[p])))*np.nan
#     Dq_P = np.empty((len(DQ_loc[p])))*np.nan
#     iN = np.where(DQ_loc[p] < 0)
#     iP = np.where(DQ_loc[p] >=0)
    
#     Dq_N[iN] = DQ_loc[p][iN]
#     Dq_P[iP] = DQ_loc[p][iP]
    
#     plt.loglog(r[p]/1000, np.abs(Dq_N), '.-r')
#     plt.loglog(r[p]/1000, np.abs(Dq_P), '.-b')
    
#     popt, pcov = curve_fit(second_order_fit, r[p]/1000, np.abs(Dq_P))
#     print(popt)

#     # plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**(2/3) + popt[1]*(bin_edges/1000)**2 - popt[2]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), '--k', linewidth=2, label='$\sim s^3$')    
#     plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**2 - popt[1]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), '--k', linewidth=2, label='$\sim s^3$')    
# plt.grid()
# plt.show()




# for i in range(num_bins):
#     mask = (r_l >= bin_edges[i]) & (r_l < bin_edges[i + 1])
#     if np.sum(mask) > 0:
#         # print(np.sum(mask))
#         dq_loc_binned[i] = np.mean(dtw[mask])

# DQ[row,:] = np.array(dq_loc_binned)         
# row = row+1







































# ### Kurtosis
# plt.figure()
# for c in [15,16,17,18]:
#     wind_speed_data = WINDS[:,c]
    
#     # Example usage:
#     q = 3  # Replace with the desired q value
#     max_tau = 288*15  # Replace with the maximum value of tau you want to calculate
#     D_2 = calculate_structure_function_q(wind_speed_data, 2, max_tau)
#     D_4 = calculate_structure_function_q(wind_speed_data, 4, max_tau)

    
    
    
#     # D_q_loc = calculate_structure_function_q_loc(wind_speed_data, q, max_tau)
    
#     # D_q contains D_q(1) to D_q(max_tau)
#     # D_q_loc contains D_q^loc(r) for each tau from 1 to max_tau
    
    
#     # plt.figure()
#     # for d in D_q_loc:
#     #     plt.plot(d,'-o')
    
#     tau = np.arange(5*60,(max_tau+1)*5*60,5*60) 
    
    
#     plt.semilogx(tau,D_4/D_2**2,'.')
    
# plt.grid()


