# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:57:02 2024

I use this code to extract the winds from MAARSY and save it in arrays to process or plot later.

@author: ma042
"""

import os
import numpy as np
import h5py as h5
from datetime import datetime

def extract_winds(ipath, year, month, component=0):
    U0 = []
    V0 = []
    U0mod = []
    heights = []
    date_list = []

    for y in year:
        for m in month:
            for d in range(1, 32):
                filename = ipath + '%04d/%04d%02d/VEL_h5/%04d%02d%02d_MST_MAARSY_Tropo_LTUT_windVEL_t15_v01c.h5' % (y, y, m, y, m, d)

                if os.path.isfile(filename):
                    with h5.File(filename, 'r') as f:
                        u = np.array(f.get('wind/u'))
                        v = np.array(f.get('wind/v'))
                        w = np.array(f.get('wind/w'))
                                                
                        if component == 0:
                            comp = u
                        elif component == 1:
                            comp = v
                        elif component == 2:
                            comp = np.sqrt(u**2 + v**2)
                        elif component == 3:
                            comp = np.sqrt(u**2 + v**2 + w**2)
                        
                        datenum = np.squeeze(f.get('info/datenums/'))
                        hours = np.squeeze(f.get('info/hour/')).astype(int)
                        mins = np.squeeze(f.get('info/min/')).astype(int)
                        secs = np.squeeze(f.get('info/sec/')).astype(int)
                        heights = np.array(f.get('info/alt'))
                        
                        # Construct the base date for each file
                        base_date = datetime(y, m, d)
                        # Construct date_list for each file
                        date_list += [base_date.replace(hour=hour, minute=minute, second=second) for hour, minute, second in zip(hours, mins, secs)]

                        U0.append(u)
                        V0.append(v)
                        U0mod.append(comp)
                        # heights.append(heights)
                
                else:
                    print(f"File not found for {y}{m:02d}{d:02d}. Skipping...")

    if U0:
        U0 = np.concatenate(U0, axis=0)
        V0 = np.concatenate(V0, axis=0)
        U0mod = np.concatenate(U0mod, axis=0)
        # heights = np.concatenate(heights, axis=0)

        return U0, V0, U0mod, np.squeeze(heights), date_list