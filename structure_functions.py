# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 17:02:38 2024

Routines to estimate structure functions from time series of winds

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


def calculate_structure_function_q_loc_binned(wind_speed_data, q, max_tau, Tau, s_min, s_max, num_bins, central_tendency='mean'):    
    n = len(wind_speed_data)
    
    bin_edges = np.logspace(np.log10(s_min), np.log10(s_max), num_bins + 1)    
    
    R_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan
    Rcount_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan
    DQ_loc_binned = np.empty((len(bin_edges),max_tau))*np.nan
    vmncount_loc_binned = []

    for ntau in range(1, max_tau + 1):
        # print(ntau)
        # v_avg = np.convolve(wind_speed_data, np.ones(ntau)/ntau, mode='valid')     
        v_avg = (1/ntau)*np.convolve(wind_speed_data, np.ones(ntau), mode='valid')     
        
        r_values = v_avg[:-1] * Tau[ntau-1] # check [:-1] or [1:] : only there to match dimensions with delta_v_q --> no afecta en nada
        # r_values = v_avg[1:] * Tau[ntau-1] 
        delta_v_q = (wind_speed_data[:n - ntau] - wind_speed_data[ntau:]) ** q
           
        vmncount_loc_binned.append(v_avg)
        
        # print(min(v_avg),max(v_avg))
        # print(min(r_values/1000),max(r_values/1000))
        # print(delta_v_q.shape)
        # print(len(r_values))
        
        # Binning in s
        for b in range(num_bins):
            ib = np.where((r_values >= bin_edges[b]) & (r_values < bin_edges[b + 1]))
            Rcount_loc_binned[b,ntau-1] = len(ib[0])

            if central_tendency=='median':        
                R_loc_binned[b,ntau-1] = np.nanmedian(r_values[ib])
                DQ_loc_binned[b,ntau-1] = np.nanmedian(delta_v_q[ib])
                
            elif central_tendency=='mean':                
                R_loc_binned[b,ntau-1] = np.nanmean(r_values[ib])
                DQ_loc_binned[b,ntau-1] = np.nanmean(delta_v_q[ib])

    return DQ_loc_binned, R_loc_binned, Rcount_loc_binned, vmncount_loc_binned





