#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:01:16 2022

@author: natsukoyamaguchi

Calculate lris spectra with donor-contribution subtracted (i.e. in-eclipse spectrum subtracted)

"""
import numpy as np
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from pathlib import Path

temp_dnr = 4250
folder_loc = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/DonorSubtractedSpectra/'
Path(folder_loc + 'outputs/').mkdir(parents=True, exist_ok=True)
save_loc = folder_loc + 'outputs/' 

## If wanting to subtract model spectra ##
    
# donor_model = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/WD/outputs/donor_' + str(temp_dnr) + '/spec_donor_unnorm.csv')
# donor_interp = interp1d(donor_model['wav'], donor_model['flux'])

## If wanting to subtract in-eclipse spectra ##

# 1. Get unnormalized spectra 
base_directory = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis'
files = sorted(listdir(base_directory + '/lris_p274')) # alphabetical order

dat_nan_b = []
dat_nan_r = []

for i in range(len(files)):
    if files[i][0:5] == 'celsb':
        dat_nan_b.append(np.loadtxt(base_directory + '/lris_p274' + '/' + files[i])) 
    elif files[i][0:5] == 'celsr':
        dat_nan_r.append(np.loadtxt(base_directory + '/lris_p274' + '/' + files[i])) 

dat_b = [] 
dat_r = []
for i in range(len(dat_nan_b)):
    dat_i = dat_nan_b[i]
    dat_b.append(dat_i[~np.isnan(dat_i).any(axis=1)]) 
for i in range(len(dat_nan_r)):
    dat_i = dat_nan_r[i]
    dat_r.append(dat_i[~np.isnan(dat_i).any(axis=1)]) 
   
# 2. Get the corresponding time (i.e. jdmid, phase) and RV shift info calculated before

r_times = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/r_times.csv')
b_times = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/b_times.csv')

r_shifts = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/outputs/shift_r.csv')
b_shifts = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/outputs/shift_b.csv')

# 4. Find spectra in each filter that is closest to eclipse 

idx = np.argwhere(np.abs(0.5-r_times['phase']).to_numpy() ==  np.min(np.abs(0.5-r_times['phase']).to_numpy()))[0][0]
idx_b = np.argwhere(np.abs(0.5-b_times['phase']).to_numpy() ==  np.min(np.abs(0.5-b_times['phase']).to_numpy()))[0][0]

# r spectra closest to eclipse, phase ~ 0.495
x_donor = dat_r[idx][:,0] - dat_r[idx][:,0] * (np.exp(r_shifts['log_shift'][29])-1)
donor_spectra = pd.DataFrame(data = {'wav' : dat_r[idx][:,0], 'wav_shift': x_donor, 'flux': dat_r[idx][:,1]})
donor_interp = interp1d(donor_spectra['wav_shift'], donor_spectra['flux'], fill_value = 'extrapolate')

# b spectra closest to eclipse, phase ~ 0.505
x_donor_b = dat_b[idx_b][:,0] - dat_b[idx_b][:,0] * (np.exp(b_shifts['log_shift'][36])-1)
donor_spectra_b = pd.DataFrame(data = {'wav' : dat_b[idx_b][:,0], 'wav_shift': x_donor_b, 'flux': dat_b[idx_b][:,1]})
donor_interp_b = interp1d(donor_spectra_b['wav_shift'], donor_spectra_b['flux'], fill_value = 'extrapolate')

for i in range(len(b_times)):
    
    b_shift = np.exp(b_shifts['log_shift'][i])-1 
    b_spectra = pd.DataFrame(data = {'wav': dat_b[i][:,0], 'wav_shift': dat_b[i][:,0] * (1 - b_shift), 'flux':dat_b[i][:,1]})
    
    # 3. Find pairs of r and b spectra that match closest in time 
    # (There is more spectra in b than r so we find a matching r for every b)

    b_jdmid = b_times['jdmid'][i]
    r_minus_b = np.abs(r_times['jdmid'] - b_jdmid)
    r_idx = r_minus_b[r_minus_b == np.min(r_minus_b)].index[0] 
    r_jdmid = r_times['jdmid'][r_idx]
    
    r_shift = np.exp(r_shifts['log_shift'][r_idx])-1 
    r_spectra = pd.DataFrame(data = {'wav': dat_r[r_idx][:,0], 'wav_shift' :dat_r[r_idx][:,0] * (1 - r_shift), 'flux':dat_r[r_idx][:,1]})

    # 4. Calculate the mean ratio between the r spectra for some wavelength range (with little atmospheric lines) and the donor (i.e. in-eclipse) spectra 
    r_slice = r_spectra[(r_spectra['wav_shift'] > 9125) & (r_spectra['wav_shift'] < 9175)] 
    donor_flux = donor_interp(r_slice['wav_shift'])
    scaling = np.mean(r_slice['flux']/donor_flux)
    
    # plt.plot(r_spectra['wav_shift'], r_spectra['flux'], label  = 'out of eclipse')
    # plt.plot(donor_spectra['wav_shift'], donor_spectra['flux'] * scaling, label  = 'in eclipse')
    # plt.legend()
    
    # 5. Subtract the flux from the scaled donor spectrum from the flux of both spectra
    r_donor_sub = pd.DataFrame(data = {'wav': r_spectra['wav'], 'wav_shift': r_spectra['wav_shift'], 
                                       'dnr_sub_flux': r_spectra['flux'] - (donor_interp(r_spectra['wav_shift']) * scaling)})
    
    b_donor_sub = pd.DataFrame(data = {'wav': b_spectra['wav'], 'wav_shift': b_spectra['wav_shift'], 
                                       'dnr_sub_flux': b_spectra['flux'] - (donor_interp_b(b_spectra['wav_shift']) * scaling)})
    
    
    # plt.plot(b_spectra['wav_shift'], b_spectra['flux'])
    # plt.plot(donor_spectra_b['wav_shift'], donor_spectra_b['flux']* scaling)
    
    b_donor_sub.to_csv(save_loc + 'b_{0:.5f}'.format(b_jdmid) + '.csv', index = False)
    r_donor_sub.to_csv(save_loc + 'r_{0:.5f}'.format(r_jdmid) + '.csv', index = False)
    
    # plt.plot(donor_model['wav'], donor_model['flux'] * scaling)


