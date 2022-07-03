#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:46:40 2022

@author: natsukoyamaguchi

Get template spectrum, with instrumental and rotational broadening applied, to use in CCF to get RV plots 

"""
from base_settings import *
from astropy.io import fits

# For rotational broadening
eps = 0.5 # linear limb-darkening coeff
rv = (2 * np.pi * 0.23 * 695700)/(P*86400)  # km/s, projected rotational velocity

#%% Normalize 

# BOSZ Kurucz Model - T_eff = 5000K, log g = 4.5, 0, 0, 0, Inst broad = 20,000
hdul = fits.open(base_directory + '/templates/amp00cp00op00t5000g45v20modrt0b20000rs.fits')

template = hdul[1].data

x_temp = template.Wavelength.byteswap().newbyteorder()
y_temp = template.SpecificIntensity.byteswap().newbyteorder()

if col == 'r':
    x_lolim = 5400
    x_uplim = 6600
elif col == 'b':
    x_lolim = 3800
    x_uplim = 5200
    
x_temp_slice = x_temp[np.logical_and(x_temp >= x_lolim, x_temp <= x_uplim)]
y_temp_slice = y_temp[np.logical_and(x_temp >= x_lolim, x_temp <= x_uplim)]
    
if norm == 'poly':
    fit_poly = np.polyfit(x_temp_slice, y_temp_slice, 2)
    polyfit_norm = fit_poly[0]*x_temp_slice**2 + fit_poly[1]*x_temp_slice + fit_poly[2]
    flux_norm_temp = y_temp_slice/polyfit_norm
elif norm == 'rolling':
    grid = np.arange(np.min(x_temp_slice), np.max(x_temp_slice), 0.14)
    interp_func = interp1d(x_temp_slice, y_temp_slice)
    rolling = pd.Series(interp_func(grid)).rolling(window, min_periods = 1).median()
    interp_func_back = interp1d(grid, rolling, fill_value="extrapolate")
    flux_norm_temp = y_temp_slice/interp_func_back(x_temp_slice)

temp_spectrum = pd.DataFrame(data = {'wavelength': x_temp_slice, 'flux': flux_norm_temp})

diff = []
for i in range(len(x_temp_slice)-1):
    diff.append(x_temp_slice[i+1]-x_temp_slice[i])
    
grid = np.mean(diff) # average wavelength separation between consecutive points used for the grid

interp_func = interp1d(temp_spectrum['wavelength'], temp_spectrum['flux'])
x = np.arange(temp_spectrum['wavelength'][0], np.max(temp_spectrum['wavelength']), grid)
y = interp_func(x)

# inst then rot 

## Instrumental Broadening (Gaussian kernel convolution) ##

from astropy.convolution import Gaussian1DKernel, convolve

# FWHM 

if norm == 'poly':
    skyline = pd.read_csv(base_directory + '/outputs' + '/' + 'skyline' + '_' + norm + '_' + col + '.csv')
elif norm == 'rolling':
    skyline = pd.read_csv(base_directory + '/outputs' + '/' + 'skyline' + '_' + norm + '_' + str(window) + '_' + col + '.csv')


FWHM = np.mean(abs(skyline['FWHM']))

sd = FWHM / (grid * 2*np.sqrt(2 * np.log(2))) # want sigma instead of FMHM
g = Gaussian1DKernel(stddev=sd)

z = convolve(y, g)

# plt.plot(x, y, 'k-', label='Before', lw = 0.8)
# plt.plot(x, z, 'b-', label='After', alpha=0.5, lw = 0.8)
# plt.legend(loc='best')
# plt.show()

temp_spectrum_inst = pd.DataFrame(data = {'log_wavelength': np.log(x), 'flux': z})

## Rotational Broadening ##

from PyAstronomy import pyasl

rflux = pyasl.rotBroad(np.exp(temp_spectrum_inst['log_wavelength']).to_numpy(), temp_spectrum_inst['flux'].to_numpy(), epsilon = eps, vsini = rv)

temp_spectrum_inst_rot = pd.DataFrame(data = {'log_wavelength': np.log(x), 'flux': rflux})

if norm == 'poly':
    temp_spectrum_inst_rot.to_csv(base_directory + '/outputs' + '/temp_inst_rot_' + norm + '_' + col + '.csv', index = False)
elif norm == 'rolling':
    temp_spectrum_inst_rot.to_csv(base_directory + '/outputs' + '/temp_inst_rot_' + norm + '_' + str(window) + '_' + col + '.csv', index = False)

#%% rot then inst 

# ## Rotational Broadening ##

# from PyAstronomy import pyasl

# rflux = pyasl.rotBroad(x, y, epsilon = eps, vsini = rv)

# temp_sectrum_rot = pd.DataFrame(data = {'wavelength': np.log(x), 'flux': rflux})

# ## Instrumental Broadening (Gaussian kernel convolution) ##

# from astropy.convolution import Gaussian1DKernel, convolve

# # FWHM 
# skyline = pd.read_csv(base_directory + 'outputs/skyline_' + norm + '_' + col + '.csv')
# FWHM = np.mean(abs(skyline['FWHM']))

# sd = FWHM / (grid * 2*np.sqrt(2 * np.log(2))) # want sigma instead of FMHM
# g = Gaussian1DKernel(stddev=sd)

# z = convolve(temp_sectrum_rot['flux'], g)

# temp_spectrum_rot_inst = pd.DataFrame(data = {'wavelength': np.log(x), 'flux': z})
# temp_spectrum_rot_inst.to_csv(base_directory + '/outputs' + '/temp_rot_inst_' + norm + '_' + col + '.csv', index = False)

# #%% 

# col = 'r'

# if col == 'r':
#     x_lolim = 5400
#     x_uplim = 6600
# elif col == 'b':
#     x_lolim = 3800
#     x_uplim = 5200
    
# temp_spectrum_inst_rot = pd.read_csv(base_directory + '/outputs' + '/temp_inst_rot_' + col + '.csv')
# temp_spectrum_rot_inst = pd.read_csv(base_directory + '/outputs' + '/temp_rot_inst_' + col + '.csv')

# rot_inst = temp_spectrum_rot_inst[(np.exp(temp_spectrum_rot_inst['wavelength']) > x_lolim + 200) & (np.exp(temp_spectrum_rot_inst['wavelength']) < x_uplim - 200)]
# inst_rot = temp_spectrum_inst_rot[(np.exp(temp_spectrum_inst_rot['wavelength']) > x_lolim + 200) & (np.exp(temp_spectrum_inst_rot['wavelength']) < x_uplim - 200)]
# plt.plot(rot_inst['wavelength'], rot_inst['flux'] - inst_rot['flux'] , lw = 0.3)

# plt.xlabel('log(wavelength)')
# plt.ylabel('flux difference')
# plt.legend()

# #%%
# plt.plot(inst_rot['wavelength'], inst_rot['flux'], lw = 0.3, label = 'Inst then Rot')
# plt.plot(rot_inst['wavelength'], rot_inst['flux'], lw = 0.3, label = 'Rot then Inst')
# plt.legend()

# plt.xlabel('log(wavelength)')
# plt.ylabel('flux')