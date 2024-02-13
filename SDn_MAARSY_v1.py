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

## self-generated codes
from extract_winds import extract_winds
from structure_functions import calculate_structure_function_q, calculate_structure_function_q_loc


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



# SAVE FILES
with open('C:/Poblet/IAP/work/N-order_sf/data/date_list_2022.txt', 'w') as f:
    f.write('\n'.join([date.strftime('%Y-%m-%d %H:%M:%S') for date in date_list]))


np.save('C:/Poblet/IAP/work/N-order_sf/data/WIND2022.npy', WINDS)
np.save('C:/Poblet/IAP/work/N-order_sf/data/height.npy', H)


# LOAD FILES
date_list = [line.strip() for line in open('C:/Poblet/IAP/work/N-order_sf/data/date_list_2022.txt')]
WINDS = np.load('C:/Poblet/IAP/work/N-order_sf/data/WIND2022.npy')
H = np.load('C:/Poblet/IAP/work/N-order_sf/data/height.npy')


### Third order SFs plots
## Altitude range:
ih0 = np.squeeze(np.where((H>=1) & (H<6.5)))
ih1 = np.squeeze(np.where((H>=6.5) & (H<11)))
ih2 = np.squeeze(np.where((H>=11) & (H<15)))
ihT = np.squeeze(np.where((H>=1) & (H<20)))

# Example usage:
q = 3  # Replace with the desired q value
max_tau = 3000  # Replace with the maximum value of tau you want to calculate
tau = np.arange(15*60,(max_tau+1)*15*60,15*60) # in seconds
i_tt = np.where((tau == 1 * 60 * 60) | (tau == 24 * 60 * 60) | (tau == 3 * 24 * 60 * 60) | (tau == 30 * 24 * 60 * 60))

for ih in [ih0, ih1, ih2]:
    fig, ax = plt.subplots()
    Dq_13 = []
    for c in ih:
        wind_speed_data = WINDS[:, c]
        D_q = calculate_structure_function_q(wind_speed_data, q, max_tau)
        # D_q = calculate_structure_function_q_loc(wind_speed_data, q, max_tau)
        d_q13 = []
        for d_q in D_q:
            d_q13.append(cube(d_q))
        Dq_13.append(np.array(d_q13))
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


