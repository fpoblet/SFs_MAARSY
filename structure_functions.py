# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 17:02:38 2024

Routines to estimate structure functions from MAARSY winds

@author: ma042
"""
import numpy as np

# routines to estimate high-order structure functions




def calculate_structure_function_q(wind_speed_data, q, max_tau):
    n = len(wind_speed_data)
    D_q = np.zeros(max_tau)

    for tau in range(1, max_tau + 1):               
        delta_v = wind_speed_data[:n - tau] - wind_speed_data[tau:]
        D_q[tau - 1] = np.nanmean(delta_v ** q)

    return D_q


def calculate_structure_function_q_loc(wind_speed_data, q, max_tau, Tau):    
    n = len(wind_speed_data)
    Dq_loc_mean = np.empty(max_tau)*np.nan
    r_loc_mean = np.empty(max_tau)*np.nan

    for ntau in range(1, max_tau + 1):
        v_avg = np.convolve(wind_speed_data, np.ones(ntau), 'valid')/ntau     
        
        r_values = v_avg * Tau[ntau-1]
        
        delta_v_q = (wind_speed_data[:n - ntau] - wind_speed_data[ntau:]) ** q

        Dq_loc_mean[ntau-1] = np.nanmean(delta_v_q)
        r_loc_mean[ntau-1] = np.nanmean(r_values)

    return Dq_loc_mean, r_loc_mean


def calculate_structure_function_q_loc_binned(wind_speed_data, q, max_tau, Tau, s_min, s_max, num_bins,median=False):    
    n = len(wind_speed_data)
    
    bin_edges = np.logspace(np.log10(s_min), np.log10(s_max), num_bins + 1)    
    
    R_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan
    Rcount_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan
    DQ_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan

    for ntau in range(1, max_tau + 1):
        # print(ntau)
        v_avg = np.convolve(wind_speed_data, np.ones(ntau)/ntau, 'valid')     
        
        r_values = v_avg[:-1] * Tau[ntau-1] # check[:-1], only there to match dimensions with delta_v_q
        delta_v_q = (wind_speed_data[:n - ntau] - wind_speed_data[ntau:]) ** q
           
        # print(min(r_values/1000),max(r_values/1000))
        # print(delta_v_q.shape)
        # print(len(r_values))
        
        for b in range(num_bins):
            ib = np.where((r_values >= bin_edges[b]) & (r_values < bin_edges[b + 1]))
            # print(len(ib[0]))
            
            if median:
                R_loc_binned[b,ntau-1] = np.nanmedian(r_values[ib])
                Rcount_loc_binned[b,ntau-1] = len(ib[0])
                DQ_loc_binned[b,ntau-1] = np.nanmedian(delta_v_q[ib])
            else:
                R_loc_binned[b,ntau-1] = np.nanmean(r_values[ib])
                Rcount_loc_binned[b,ntau-1] = len(ib[0])
                DQ_loc_binned[b,ntau-1] = np.nanmean(delta_v_q[ib])

    return DQ_loc_binned, R_loc_binned, Rcount_loc_binned



# def calculate_structure_function_q_loc(wind_speed_data, q, max_tau):
#     n = len(wind_speed_data)
#     Dq_loc_mean = []
#     Dq_loc_median = []
#     r_loc_mean = []
#     r_loc_median = []

#     for tau in range(1, max_tau + 1):
       
#         r_values = []
#         delta_v_q = []

#         for t in range(n - tau):
#             v_avg = np.mean(wind_speed_data[t:t + tau])
#             r = v_avg * tau
#             r_values.append(r)
#             delta_v_q.append((wind_speed_data[t] - wind_speed_data[t + tau]) ** q)

#         r_values_mean = np.nanmean(r_values)
#         r_values_median = np.nanmedian(r_values)
#         r_values_std = np.nanstd(r_values)
#         delta_v_q_mean = np.nanmean(delta_v_q)
#         delta_v_q_median = np.nanmedian(delta_v_q)

#         # # Calculate D_q^loc(r) by averaging in fifty bins every decade
#         # r_min = r_values.min()
#         # r_max = r_values.max()
#         # # num_bins = 50
#         # # bin_edges = np.logspace(np.log10(r_min), np.log10(r_max), num_bins + 1)
#         # D_q_loc_tau = []

#         # for i in range(num_bins):
#         #     mask = (r_values >= bin_edges[i]) & (r_values < bin_edges[i + 1])
#         #     if np.sum(mask) > 0:
#         #         D_q_loc_tau.append(np.mean(delta_v_q[mask]))

#         Dq_loc_mean.append(delta_v_q_mean)
#         Dq_loc_median.append(delta_v_q_median)
#         r_loc_mean.append(r_loc_mean)
#         r_loc_median.append(r_loc_median)

#     return np.array(Dq_loc_mean), np.array(Dq_loc_mean), np.array(r_loc_mean), np.array(r_loc_median)
