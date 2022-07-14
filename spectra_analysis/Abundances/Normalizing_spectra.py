#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 12:01:25 2022

@author: natsukoyamaguchi
"""

from Which_spectra import *

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances'

Path(folder_dir + '/spectra_broadened').mkdir(parents=True, exist_ok=True)

if norm == 'poly':
    file_dir = folder_dir + '/spectra_broadened' + '/' + norm + '/' + col
    save_dir = folder_dir + '/spectra_norm' + '/' + norm + '/' + col
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    files = sorted(listdir(file_dir))
elif norm == 'rolling':
    file_dir = folder_dir + '/spectra_broadened' + '/' + norm + '/' + str(window) +  '/' + col
    save_dir = folder_dir + '/spectra_norm' + '/' + norm + '/' + str(window) + '/' + col
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    files = sorted(listdir(file_dir)) 

for i in range(len(files)):
    dat = pd.read_csv(file_dir + '/' + files[i])
    if norm == 'poly':
        z = np.polyfit(dat['wavelength'], dat['flux'], 2)
        fit = z[0] * dat['wavelength']**2 + z[1] * dat['wavelength'] + z[2]
        # plt.plot(dat['wavelength'], dat['flux'])
        # plt.plot(dat['wavelength'], fit)
        spectrum_norm = pd.DataFrame(data = {'wavelength': dat['wavelength'], 'flux' : dat['flux']/fit})
        spectrum_norm.to_csv(save_dir + '/' + files[1:len(files)][i] + '.csv', index = False)
    elif norm == 'rolling':
        grid = np.arange(np.min(dat['wavelength']), np.max(dat['wavelength']), 0.002)
        interp_func = interp1d(dat['wavelength'], dat['flux'], fill_value="extrapolate")
        rolling = pd.Series(interp_func(grid)).rolling(int(window/0.002), min_periods = 1, center = True).median()
        interp_func_back = interp1d(grid, rolling, fill_value="extrapolate")
        cont = interp_func_back(dat['wavelength'])
        # plt.plot(dat['wavelength'], dat['flux'])
        # plt.plot(dat['wavelength'], cont)
        spectrum_norm = pd.DataFrame(data = {'wavelength': dat['wavelength'], 'flux' : dat['flux']/cont}) # convert wavelength into AA
        spectrum_norm.to_csv(save_dir + '/' + files[i] + '.csv', index = False)
    