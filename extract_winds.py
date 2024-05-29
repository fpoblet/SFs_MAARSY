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
from datetime import datetime, timedelta

import matplotlib.pyplot as plt


# FUNCTIONS
# Function to find positions of common elements in seq2
def find_common_elements_positions(seq1, seq2):
    positions = [seq2.index(dt) for dt in seq1 if dt in seq2]
    return positions

def extract_winds(ipath, year, month, time_sequence, component=2):
    # U0 = []
    # V0 = []
    # U0mod = []
    # heights = []

    date_list = []

    U0 = np.empty((len(time_sequence),77))*np.nan # 77 number of height levels
    V0 = np.empty((len(time_sequence),77))*np.nan
    U0mod = np.empty((len(time_sequence),77))*np.nan
    heights = np.empty((len(time_sequence),77))*np.nan

    for y in year:
        for m in month:
            for d in range(1, 32):
                filename = ipath + '%04d/%04d%02d/VEL_h5/%04d%02d%02d_MST_MAARSY_Tropo_LTUT_windVEL_t15_v01c.h5' % (y, y, m, y, m, d)

                if os.path.isfile(filename):
                    with h5.File(filename, 'r') as f:
                        u = np.array(f.get('wind/u'))
                        v = np.array(f.get('wind/v'))
                        w = np.array(f.get('wind/w'))
                                                
                        # if component == 0:
                        #     comp = u
                        # elif component == 1:
                        #     comp = v
                        # elif component == 2:
                        #     comp = np.sqrt(u**2 + v**2)
                        # elif component == 3:
                        #     comp = np.sqrt(u**2 + v**2 + w**2)
                        
                        datenum = np.squeeze(f.get('info/datenums/'))
                        hours = np.squeeze(f.get('info/hour/')).astype(int)
                        mins = np.squeeze(f.get('info/min/')).astype(int)
                        secs = np.squeeze(f.get('info/sec/')).astype(int)
                        
                        heights = np.array(f.get('info/alt'))
                        
                        # Construct the base date for each file
                        base_date = datetime(y, m, d)

                        # # Construct date_list for each file
                        # date_list += [base_date.replace(hour=hour, minute=minute, second=second) for hour, minute, second in zip(hours, mins, secs)]
                        # print(date_list)

                        date_list = [base_date.replace(hour=hour, minute=minute, second=second) for hour, minute, second in zip(hours, mins, secs)]
                        # print(date_list)

                        # Find common elements and their positions in time_sequence
                        positions = find_common_elements_positions(date_list, time_sequence)
                        # print(f"{y}{m:02d}{d:02d}")
                        # print('date_list', len(date_list))
                        # print('positions', len(positions))
                        # print('u len',len(u))

                        U0[positions,:] = u
                        V0[positions,:] = v
                        U0mod[positions,:] = np.sqrt(u**2 + v**2)
                        
                        # U0.append(u)
                        # V0.append(v)
                        # U0mod.append(comp)
                        # # heights.append(heights)
                
                else:
                    print(f"File not found for {y}{m:02d}{d:02d}. Skipping...")

    # if U0:
    #     U0 = np.concatenate(U0, axis=0)
    #     V0 = np.concatenate(V0, axis=0)
    #     U0mod = np.concatenate(U0mod, axis=0)
    #     # heights = np.concatenate(heights, axis=0)

        # return U0, V0, U0mod, np.squeeze(heights), date_list
    return U0, V0, U0mod, np.squeeze(heights)
    
# DEFINE DATES
year = [2019, 2020, 2021, 2022]
months = [1,2,3,4,5,6,7,8,9,10,11,12]
# dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]

# year = [2021, 2022]
# months = [1,2,3,4,5,6,7,8,9,10,11,12]
# # dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]

# # SUMMERS
# year = [2019, 2020, 2021, 2022]
# months = [6,7,8]

# # WINTERS
# year = [2019, 2020, 2021, 2022]
# months = [11,12,1]

# for y in [2019, 2020, 2021,2022]:
#     year = [y]
#     months = [1,2,3,4,5,6,7,8,9,10,11,12]
#     dayssn = [31,28,31,30,31,30,31,31,30,31,30,31]


# Generate a 1-min CONTINOUS time sequence to be filled with the data, using "extract_winds" function
## Define the start and end dates
start_date = datetime(year[0], 1, 1)
end_date = datetime(year[-1], 12, 31, 23, 45)

# Generate the time sequence with 1-minute resolution
time_sequence = []
current_time = start_date
delta = timedelta(minutes=1)

while current_time <= end_date:
    time_sequence.append(current_time)
    current_time += delta

# # Display the first few elements of the sequence to verify
# for timestamp in time_sequence[:10]:  # Display the first 10 timestamps
#     print(timestamp)


####################################################################################################
# CALL FUNCTION TO EXTRACT DATA
path = 'C:/Poblet/IAP/work/N-order_sf/MAARSY/'

# U0, V0, WINDS, H, date_list =  extract_winds(ipath = path,
#                                               year = year,
#                                               month = months,
#                                               time_sequence = time_sequence,
#                                               component = 2) # 0: zonal, 1:meridional, 2:horizontal absolute value, 3: full amplitude

U0, V0, U0mod, H =  extract_winds(ipath = path,
                                  year = year,
                                  month = months,
                                  time_sequence = time_sequence,
                                  component = 2) # 0: zonal, 1:meridional, 2:horizontal absolute value, 3: full amplitude

####################################################################################################
# REDUCE THE RESOLUTION TO 15-min to work with smaller arrays and lists
# Generate the new time sequence with 15-minute resolution
time_sequence_15min = []
current_time = start_date
delta = timedelta(minutes=15)

while current_time <= end_date:
    time_sequence_15min.append(current_time)
    current_time += delta

# Initialize an empty array to store the averaged values
U0mod_15min = np.empty((len(time_sequence_15min), U0mod.shape[1]))*np.nan

# Iterate over the new time sequence
for i, new_time in enumerate(time_sequence_15min):
    # Find the corresponding indices in the original time sequence for the 15-minute window
    start_index = int((new_time - start_date).total_seconds() / 60)
    end_index = min(start_index + 15, len(time_sequence))  # Ensure not to go out of bounds

    # Calculate the average of the corresponding values in the original array
    U0mod_15min[i] = np.nanmean(U0mod[start_index:end_index], axis=0)







# # Create a plot
# plt.figure(figsize=(12, 6))
# # Use pcolormesh or imshow to create the colormap
# plt.pcolormesh(U0mod_15min.T, cmap='viridis', shading='auto')
# # Set x-ticks to match the time sequence
# plt.xticks(ticks=np.arange(len(time_sequence_15min)), labels=[dt.strftime('%Y-%m-%d %H:%M') for dt in time_sequence_15min], rotation=45, ha='right')

# # Set labels and title
# plt.xlabel('Time')
# plt.ylabel('Data Index')
# plt.title('Colormap of Averaged Data Array')

# # Show color bar
# plt.colorbar(label='Value')

# # Tight layout for better spacing
# plt.tight_layout()

# # Show plot
# plt.show()




plt.figure(figsize=(10, 6))
plt.pcolormesh(U0mod_15min.T, shading='auto')


  
####################################################################################################
# SAVE FILES
    
# Save the new time sequence to a txt file
with open('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/date_list_15min_2019-2022.txt', 'w') as file:
    for time in time_sequence_15min:
        file.write(f"{time}\n") 

np.save('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/U0mod_15min_2019-2022.npy', U0mod_15min)
np.save('C:/Poblet/IAP/work/N-order_sf/outputs_extract_winds/height.npy', H)
    


    
    
    
    
    