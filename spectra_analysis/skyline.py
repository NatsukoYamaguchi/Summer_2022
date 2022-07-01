#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:20:55 2022

@author: natsukoyamaguchi

Fit gaussian to the OI line in the sky spectra to obtain a FWHM 

"""

from base_settings import *

# Get the lris spectrum 

times = pd.read_csv(base_directory + '/spectra_norm/' + col + '_' + 'times' + '.csv')

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
            dat_i = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + files[i])
        elif norm == 'rolling':
            dat_i = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + str(window) + '/' + files[i])
        dat.append(dat_i)
        
# Fit gaussian to sky line to see if the spectra uses air or vacuum wavelegnths 

centre = 5576 # ~ OI line (air: 5575.79 AA, vacuum: 5577.34 AA)

def Gaussian(z, a, b, c, d):
    return a * np.exp(-((z-b)**2)/(2*c**2)) + d

params_FWHM_list = []

for n in range(len(dat)):
    # isolate +- 20 AA around the line
    sky_line = dat[n][(dat[n]['wavelength'] >= centre - 20) & (dat[n]['wavelength'] <= centre + 20)]
    params, params_covariance = optimize.curve_fit(Gaussian, sky_line['wavelength'], sky_line['sky_flux_norm'], p0=[12 , centre ,10, 1])
    params_FWHM = np.insert(params, 4, 2 * params[2] * np.sqrt(2 * np.log(2)))
    params_FWHM_list.append(params_FWHM)

params_FWHM_df = pd.DataFrame(data = params_FWHM_list, columns = ['a', 'b', 'c', 'd', 'FWHM'])

print('average wavelength at line centre: ', np.mean(params_FWHM_df['b']))
print('average FWHM: ', np.mean(abs(params_FWHM_df['FWHM'])))

if norm == 'poly':
    params_FWHM_df.to_csv(base_directory + '/outputs' + '/' + 'skyline' + '_' + norm + '_' + col + '.csv', index = False)
elif norm == 'rolling':
    params_FWHM_df.to_csv(base_directory + '/outputs' + '/' + 'skyline' + '_' + norm + '_' + str(window) + '_' + col + '.csv', index = False)

#%% Plotting spectrum
    
# n = 60

# plt.plot(dat[n]['wavelength'], dat[n]['flux_norm'], lw = 0.3, c = '#24BA03', label = "target spectra")
# plt.plot(dat[n]['wavelength'], dat[n]['sky_flux_norm'], lw = 0.3, c = '#029AD8', label = "sky spectra")

# plt.vlines([5575.79, 5577.34], -2, 15, ls = '--', lw = 0.3)
# plt.xlabel('Wavelength')
# plt.ylabel('Flux')
# plt.legend()

# #%% Look at OI line and gaussian fit

# sky_line = dat[n][(dat[n]['wavelength'] >= centre - 20) & (dat[n]['wavelength'] <= centre + 20)]
# plt.scatter(sky_line['wavelength'], sky_line['sky_flux_norm'], s = 3)
# x = np.linspace(centre - 20, centre + 20, 200)
# plt.plot(x, Gaussian(x, params_FWHM_df['a'][n], params_FWHM_df['b'][n], params_FWHM_df['c'][n], params_FWHM_df['d'][n]))

# #%% See if line centre and FWHM changes over observations

# plt.plot(params_FWHM_df.index, params_FWHM_df['b'])
# plt.plot(jdmid, params_FWHM_df['b'])
# plt.plot(jdmid, abs(params_FWHM_df['FWHM']))
# plt.xlabel('JD_MID')
