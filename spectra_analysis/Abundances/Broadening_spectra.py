#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 08:29:29 2022

@author: natsukoyamaguchi
"""
from Which_spectra import *

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances'
files = sorted(listdir(folder_dir + '/spectra'))
files = files[1:len(files)] # alphabetical order

Path(folder_dir + '/spectra_broadened/' + norm).mkdir(parents=True, exist_ok=True)   

if col == 'r':
    x_lolim = 5400 # wavelength range used to match with template 
    x_uplim = 10000 # usually 660 but I want to see Na line at 890
elif col == 'b':
    x_lolim = 3800
    x_uplim = 5200
    
if norm == 'poly':
    save_dir = folder_dir + '/spectra_broadened' + '/' + norm + '/' + col
    Path(save_dir).mkdir(parents=True, exist_ok=True)
elif norm == 'rolling':
    save_dir = folder_dir + '/spectra_broadened' + '/' + norm + '/' + str(window) +  '/' + col
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    

skyline = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/outputs/skyline_rolling_300_' + col + '.csv')
FWHM = np.mean(abs(skyline['FWHM']))

# For rotational broadening
P = 0.11414278 #days 
eps = 0.5 # linear limb-darkening coeff
rv = (2 * np.pi * 0.23 * 695700)/(P*86400)  # km/s, projected rotational velocity

for i in range(len(files)):
    dat = np.loadtxt(folder_dir + '/spectra' + '/' + files[i])
    spectrum_full = pd.DataFrame(data = {'wavelength': dat[:,0] * 10, 'flux' : dat[:,1]})
    spectrum_df = spectrum_full[(spectrum_full['wavelength'] > x_lolim) & (spectrum_full['wavelength'] < x_uplim)]
    grid = np.mean(np.diff(spectrum_df['wavelength']))
    interp_func = interp1d(spectrum_df['wavelength'], spectrum_df['flux'], fill_value = 'extrapolate')
    x = np.arange(np.min(spectrum_df['wavelength']), np.max(spectrum_df['wavelength']), grid)
    y = interp_func(x)
    sd = FWHM / (grid * 2*np.sqrt(2 * np.log(2))) # want sigma instead of FMHM
    g = Gaussian1DKernel(stddev=sd)
    z = convolve(y, g)
    temp_spectrum_inst = pd.DataFrame(data = {'wavelength': x, 'flux': z})
    rflux = pyasl.rotBroad(temp_spectrum_inst['wavelength'].to_numpy(), temp_spectrum_inst['flux'].to_numpy(), epsilon = eps, vsini = rv)
    temp_spectrum_inst_rot = pd.DataFrame(data = {'wavelength': x, 'flux': rflux})
    temp_spectrum_inst_rot.to_csv(save_dir + '/' + files[i] + '.csv', index = False)
    # plt.plot(spectrum_df['wavelength'], spectrum_df['flux'])
    # plt.plot(temp_spectrum_inst_rot['wavelength'], temp_spectrum_inst_rot['flux'])
