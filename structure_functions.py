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
    Dq_loc_median = np.empty(max_tau)*np.nan
    r_loc_mean = np.empty(max_tau)*np.nan
    r_loc_median = np.empty(max_tau)*np.nan

    for tau in range(1, max_tau + 1):
        v_avg = np.convolve(wind_speed_data, np.ones(tau)/tau, 'valid')     
        
        r_values = v_avg * Tau[tau-1]
        
        delta_v_q = (wind_speed_data[:n - tau] - wind_speed_data[tau:]) ** q

        Dq_loc_mean[tau-1] = np.nanmean(delta_v_q)
        Dq_loc_median[tau-1] = np.nanmedian(delta_v_q)
        r_loc_mean[tau-1] = np.nanmean(r_values)
        r_loc_median[tau-1] = np.nanmedian(r_values)

    return Dq_loc_mean, Dq_loc_median, r_loc_mean, r_loc_median


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
