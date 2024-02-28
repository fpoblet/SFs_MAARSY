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
from structure_functions import calculate_structure_function_q, calculate_structure_function_q_loc


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


# year = [2019, 2020, 2021,2022]
# months = [1,2,3,4,5,6,7,8,9,10,11,12]
# dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]

year = [2022]
months = [1,2,3,4,5,6,7,8,9,10,11,12]
dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]


path = 'C:/Poblet/IAP/work/N-order_sf/MAARSY/'

U0, V0, WINDS, H, date_list =  extract_winds(ipath = path,
                                             year = year,
                                             month=months,
                                             daysn=dayssn)


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


#####################################################################################################
# SET PARAMETERS:
q = 3  # Replace with the desired q value
max_tau = 3000  # Replace with the maximum value of tau you want to calculate
tau = np.arange(15*60,(max_tau+1)*15*60,15*60) # in seconds
i_tt = np.where((tau == 1 * 60 * 60) | (tau == 24 * 60 * 60) | (tau == 3 * 24 * 60 * 60) | (tau == 30 * 24 * 60 * 60))

r_min = 10 * 1000 # minimum value in meters. It is constrained by MAARSY measurement area.
r_max = 10000 * 1000
num_bins = 500
bin_edges = np.logspace(np.log10(r_min), np.log10(r_max), num_bins + 1) # For local TH

#####################################################################################################

#####################################################################################################
### Third order SFs plots
#####################################################################################################
# TEMPORAL SFs
## Altitude range: These ranges allow to clearly distinguish different behaviors
## in termporal third-order SFs.
ih0 = np.squeeze(np.where((H>=6) & (H<6.5)))
ih1 = np.squeeze(np.where((H>=6.5) & (H<11)))
ih2 = np.squeeze(np.where((H>=11) & (H<15)))
ihT = np.squeeze(np.where((H>=1) & (H<20)))


# for ih in [ih0, ih1, ih2]:
#     print(ih)
#     fig, ax = plt.subplots()
#     Dq_13 = []
#     DQ_loc = []
#     r = []
#     for c in ih:
#         print(c)
#         wind_speed_data = WINDS[:, c]
#         D_q = calculate_structure_function_q(wind_speed_data, q, max_tau)
       
#         d_q13 = []
#         for d_q in D_q:
#             d_q13.append(cube(d_q))
#         Dq_13.append(np.array(d_q13))

#         DQ_loc.append(Dq_13)

#         # plt.semilogx(tau, d_q13, '.', label='', alpha=0.1)
    
#     # Calculate mean and standard deviation
#     mean_Dq_13 = np.nanmean(np.vstack(Dq_13), axis=0)
#     median_Dq_13 = np.nanmedian(np.vstack(Dq_13), axis=0)
#     std_Dq_13 = np.nanstd(np.vstack(Dq_13), axis=0)
    
#     # Plot mean and fill_between
#     plt.semilogx(tau, mean_Dq_13, '.k', label='Mean')
#     plt.semilogx(tau, median_Dq_13, 'xr', label='Median')
#     plt.fill_between(tau, mean_Dq_13 + std_Dq_13, mean_Dq_13 - std_Dq_13, color='gray', alpha=0.3, label='± 1 std')
    
#     plt.ylim([-10, 10])
#     plt.legend()
#     plt.grid()
    
#     plt.xlabel(r'$\tau$ (s)')
#     plt.ylabel(r'$D_3^{1/3}$ (m/s)')
#     plt.title('Heights: %1.1f-%1.1f km'%(min(H[ih]),max(H[ih])))


#     # Add vertical lines and text annotations
#     for idx in i_tt[0]:
#         plt.axvline(x=tau[idx], color='k', linestyle='--', linewidth=2.0)
#         plt.text(tau[idx], -9.7, ['1 h', '1 d', '3 d', '30 d'][i_tt[0].tolist().index(idx)], color='k', rotation=90,
#                   verticalalignment='bottom', horizontalalignment='right')

#     plt.show()
#####################################################################################################


#####################################################################################################
# SPATIAL THIRD ORDER SFs
ihs =  np.squeeze(np.where((H>=8) & (H<12)))
# ihs =  np.squeeze(np.where((H>=1) & (H<6.5)))


DQ = np.empty((len(ihs),len(bin_edges)))*np.nan
row = 0
for ih in [ihs]:    
    print(ih)

    dq_loc_binned = np.empty(len(bin_edges))*np.nan
    DQ_loc = []

    r = []
    for c in ih:
        print(c)
        wind_speed_data = WINDS[:, c]
        Dq_loc_mean, Dq_loc_median, r_loc_mean, r_loc_median = calculate_structure_function_q_loc(wind_speed_data, 
                                                                                                  q, 
                                                                                                  max_tau,
                                                                                                  tau)
        
        # With which array I am working
        dtw = Dq_loc_mean
        r_l = r_loc_mean

        DQ_loc.append(dtw)
        r.append(r_l)

## BINNED QUANTITIES



fig, ax = plt.subplots()
plt.loglog(bin_edges/1000, 0.5*(bin_edges/1000), '--k', linewidth=2, label='$\sim s^3$')
# plt.loglog(bin_edges/1000, 0.000000005*(bin_edges/1000)**3, '--k', linewidth=2, label='$\sim s^3$')

for p in range(len(r)):
    
    # Positive - Negative discrimination
    Dq_N = np.empty((len(DQ_loc[p])))*np.nan
    Dq_P = np.empty((len(DQ_loc[p])))*np.nan
    iN = np.where(DQ_loc[p] < 0)
    iP = np.where(DQ_loc[p] >=0)
    
    Dq_N[iN] = DQ_loc[p][iN]
    Dq_P[iP] = DQ_loc[p][iP]
    
    plt.loglog(r[p]/1000, np.abs(Dq_N), '.-r')
    plt.loglog(r[p]/1000, np.abs(Dq_P), '.-b')
plt.grid()
plt.show()

################################################################################
# SPATIAL SECOND ORDER SFs
ihs =  np.squeeze(np.where((H>=9) & (H<11)))
# ihs =  np.squeeze(np.where((H>=1) & (H<6.5)))
q = 2

DQ = np.empty((len(ihs),len(bin_edges)))*np.nan
row = 0
for ih in [ihs]:    
    print(ih)

    dq_loc_binned = np.empty(len(bin_edges))*np.nan
    DQ_loc = []

    r = []
    for c in ih:
        print(c)
        wind_speed_data = WINDS[:, c]
        Dq_loc_mean, Dq_loc_median, r_loc_mean, r_loc_median = calculate_structure_function_q_loc(wind_speed_data, 
                                                                                                  q, 
                                                                                                  max_tau,
                                                                                                  tau)
        
        # With which array I am working
        dtw = Dq_loc_mean
        r_l = r_loc_mean

        DQ_loc.append(dtw)
        r.append(r_l)

## BINNED QUANTITIES



fig, ax = plt.subplots()
for p in range(len(r)):
    
    # Positive - Negative discrimination
    Dq_N = np.empty((len(DQ_loc[p])))*np.nan
    Dq_P = np.empty((len(DQ_loc[p])))*np.nan
    iN = np.where(DQ_loc[p] < 0)
    iP = np.where(DQ_loc[p] >=0)
    
    Dq_N[iN] = DQ_loc[p][iN]
    Dq_P[iP] = DQ_loc[p][iP]
    
    plt.loglog(r[p]/1000, np.abs(Dq_N), '.-r')
    plt.loglog(r[p]/1000, np.abs(Dq_P), '.-b')
    
    popt, pcov = curve_fit(second_order_fit, r[p]/1000, np.abs(Dq_P))
    print(popt)

    # plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**(2/3) + popt[1]*(bin_edges/1000)**2 - popt[2]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), '--k', linewidth=2, label='$\sim s^3$')    
    plt.loglog(bin_edges/1000, popt[0]*(bin_edges/1000)**2 - popt[1]*(bin_edges/1000)**2*(np.log(bin_edges/1000)), '--k', linewidth=2, label='$\sim s^3$')    
plt.grid()
plt.show()




# for i in range(num_bins):
#     mask = (r_l >= bin_edges[i]) & (r_l < bin_edges[i + 1])
#     if np.sum(mask) > 0:
#         # print(np.sum(mask))
#         dq_loc_binned[i] = np.mean(dtw[mask])

# DQ[row,:] = np.array(dq_loc_binned)         
# row = row+1
























## Checking winds
date_nums = date2num(date_list)  # Convert to matplotlib date format

# Plot
plt.figure(figsize=(10, 6))
plt.pcolormesh(date_nums, H, WINDS.T, shading='auto')
plt.colorbar(label='Wind Velocity Modulus (m/s)')
plt.xlabel('Date')
plt.ylabel('Height (km)')
plt.title('Wind Velocity Modulus')
plt.gca().xaxis_date()  # Set the x-axis as dates
# plt.gca().invert_yaxis()  # Invert y-axis to have higher altitude at the top
plt.tight_layout()
plt.show()


plt.figure()
plt.plot(date_nums,WINDS,'.')
plt.gca().xaxis_date()


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


