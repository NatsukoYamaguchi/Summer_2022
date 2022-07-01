#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 18:44:46 2022

@author: natsukoyamaguchi

Get river plot of the spectra  

"""
from base_settings import *

#%% Settings

# wavelength range to plot 
# r: 5500 - 7000 (H alpha: 6500 - 6600, Na: 5850 - 5950)
# b: 4000 - 5500
if col == 'r':
    lo = 5500 
    hi = 7000
elif col == 'b':
    lo = 4000 
    hi = 5500
    
# vmin and vmax (imshow data range)
vmin= 0.2 # 0.5
vmax= 1.3 # 1.5

#%% Get the lris spectrum 

save_dir = '/Users/natsukoyamaguchi/Desktop/Summer 2022/lris_p274_analysis/outputs'

times = pd.read_csv(base_directory + '/spectra_norm' + '/' + col + '_' + 'times' + '.csv')

phase_list = times['phase']
jdmid = times['jdmid']
ra = times['ra']
dec = times['dec']

if norm == 'poly':
    files = sorted(listdir(base_directory + '/spectra_norm/' + norm)) # alphabetical order
elif norm == 'rolling':
    files = sorted(listdir(base_directory + '/spectra_norm/' + norm + '/' + str(window)))
    
dat = []
for i in range(len(files)):
    if files[i][0] == col:
        if norm == 'poly':
            dat_i = data = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + files[i])
        elif norm == 'rolling':
            dat_i = data = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + str(window) + '/' + files[i])
        dat.append(dat_i)

#%% 
    
diff = []
for n in range(len(dat)):
    x = dat[n]['wavelength']
    for i in range(len(x)-1):
        diff.append(x[i+1] - x[i])

# choose a range to plot 
if col == 'r':
    grid = np.arange(lo, hi, np.mean(diff))
    if norm == 'poly':
        title = 'red, ' + norm + ', ' + str(lo) + ' - ' + str(hi) + ' AA'
    elif norm == 'rolling':
        title = 'red, ' + norm + ' ' + str(window) + ', ' + str(lo) + ' - ' + str(hi) + ' AA'
elif col == 'b':
    grid = np.arange(lo, hi, np.mean(diff))
    if norm == 'poly':
        title = 'blue, ' + norm + ', ' + str(lo) + ' - ' + str(hi) + ' AA'
    elif norm == 'rolling':
        title = 'blue, ' + norm + ' ' + str(window) + ', ' + str(lo) + ' - ' + str(hi) + ' AA'

flux_interp = []
flux_min = [] 
flux_max = [] # helps in choosing vmin and vmax 
for n in range(len(dat)):
    x = dat[n]['wavelength']
    y = dat[n]['flux_norm']
    f = interp1d(x, y)
    flux_interp.append(f(grid))
    flux_min.append(np.min(f(grid)))
    flux_max.append(np.max(f(grid)))

#%% River plot using imshow

# Plot multiple periods 

multiple_periods = np.vstack((flux_interp, flux_interp, flux_interp))

plt.imshow(multiple_periods, extent = [np.min(grid), np.max(grid), np.min(phase_list), 3*np.max(phase_list)], cmap = 'gray', aspect = 100, vmin= vmin, vmax= vmax)

plt.xlabel('Wavelength', fontsize = 9)
plt.ylabel('Orbital Phase', fontsize = 9)
plt.title(title)
if norm == 'poly':
    plt.savefig(base_directory + '/plots/' + 'river' + '_' + col + '_' + norm + '_' + str(lo) + '_' + str(hi) + '.png')
if norm == 'rolling':
    plt.savefig(base_directory + '/plots/' + 'river' + '_' + col + '_' + norm + '_' + str(window) + '_' + str(lo) + '_' + str(hi) + '.png')