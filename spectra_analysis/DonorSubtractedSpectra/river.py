#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:56:07 2022

@author: natsukoyamaguchi

Make river plot of the donor-subtracted spectra 

"""
import numpy as np
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from pathlib import Path

col = 'r'

# r: 5500 - 7000 (H alpha: 6500 - 6600, Na: 5850 - 5950)
if col == 'r':
    lo = 5850 
    hi = 5950
elif col == 'b':
    lo = 4000 #4000 
    hi = 5500 #5500

folder_loc = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/DonorSubtractedSpectra/plots/'
Path(folder_loc + 'river').mkdir(parents=True, exist_ok=True)
save_loc = folder_loc + 'river/' 

# Get spectra

times = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/' + col + '_' + 'times' + '.csv')
phase_list = times['phase']
jdmid = times['jdmid']

base_directory = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/DonorSubtractedSpectra/outputs'
files = sorted(listdir(base_directory)) # alphabetical order

dat = []
for i in range(len(files)):
    if files[i][0] == col:
        dat_i = pd.read_csv(base_directory + '/' + files[i])
        dat.append(dat_i)

 
diff = []
for n in range(len(dat)):
    x = dat[n]['wav']
    for i in range(len(x)-1):
        diff.append(x[i+1] - x[i])

# choose a range to plot 
if col == 'r':
    grid = np.arange(lo, hi, np.mean(diff))
elif col == 'b':
    grid = np.arange(lo, hi, np.mean(diff))

flux_interp = []
flux_min = [] 
flux_max = [] # helps in choosing vmin and vmax 
for n in range(len(dat)):
    x = dat[n]['wav']
    y = dat[n]['dnr_sub_flux']
    f = interp1d(x, y)
    flux_interp.append(f(grid))
    flux_min.append(np.min(f(grid)))
    flux_max.append(np.max(f(grid)))

print(np.min(flux_min))
print(np.max(flux_max))

# vmin and vmax (imshow data range)
vmin= np.min(flux_min) 
vmax= np.max(flux_max) 

# Plot multiple periods 

multiple_periods = np.vstack((flux_interp, flux_interp, flux_interp))

plt.imshow(multiple_periods, extent = [np.min(grid), np.max(grid), np.min(phase_list), 3*np.max(phase_list)], cmap = 'gray', aspect = 100, vmin= vmin, vmax= vmax)

plt.xlabel('Wavelength', fontsize = 9)
plt.ylabel('Orbital Phase', fontsize = 9)

# plt.savefig(save_loc + col + '_' + str(lo) + '_' + str(hi) + '.png', dpi = 300)