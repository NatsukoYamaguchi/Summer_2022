#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 09:06:00 2022

@author: natsukoyamaguchi
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from pathlib import Path
from Spectra_conversions import get_koester_spectrum_cgs_units, kurucz_spectrum_at_a_given_distance

temp_dnr = 4750
folder_loc = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/WD/'
Path(folder_loc + 'outputs/' + 'donor_' + str(temp_dnr)).mkdir(parents=True, exist_ok=True)
save_loc = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/WD/outputs/' + 'donor_' + str(temp_dnr) + '/'

#%%

# Birth of ELMs paper
M_WD = 0.62 
R_protoWD = 0.23

# Gaia eDR3 source ID: 861540207303947776 
parallax_mas = 2.935117453431959 
dist_pc = 1000/parallax_mas

# Koester WD model - T_eff = 20,000K, log g = 8.0
wl_WD, flux_WD = get_koester_spectrum_cgs_units(path = folder_loc + 'models/koester_wd_20kK_logg8.dat.txt',
                                                d_pc = dist_pc, 
                                                mass_assume=M_WD)

# BOSZ Kurucz Model - T_eff = temp_dnr, log g = 4.5, 0, 0, 0, Inst broad = 20,000
wl_donor, flux_donor, cont_donor = kurucz_spectrum_at_a_given_distance(path = folder_loc + 'models/amp00cp00op00t' + str(temp_dnr) + 'g45v20modrt0b20000rs.fits',
                                                                       radius_Rsun=R_protoWD, 
                                                                       parallax=parallax_mas)

# full spectra with flux in physical units (erg s-1 cm-2 AA-1)
WD = pd.DataFrame(data = {'wav': wl_WD, 'flux': flux_WD})
donor = pd.DataFrame(data = {'wav': wl_donor, 'flux': flux_donor})

# isolate where the two spectra overlap and sum contribution 
WD_overlap = WD[WD['wav'] > np.min(donor['wav']) - 1] # - 1 avoid out of bounds interpolation
donor_overlap = donor[donor['wav'] < np.max(WD['wav'])]
interp_func = interp1d(WD_overlap['wav'], WD_overlap['flux'])
summed_overlap = pd.DataFrame(data = {'wav' : donor_overlap['wav'], 'flux': interp_func(donor_overlap['wav']) + donor_overlap['flux']})

summed_spectra = pd.concat([WD[WD['wav'] <= np.min(donor['wav'])], 
                            summed_overlap, 
                            donor[donor['wav'] >= np.max(WD['wav'])]])


summed_spectra = summed_spectra[(summed_spectra['wav'] < 20000)]
donor_spectra = donor[(donor['wav'] < 20000)]

#%% 

# Interpolate to a regular gird 

grid_size = 1
const_grid = np.arange(np.min(summed_spectra['wav']), np.max(summed_spectra['wav']), grid_size)
interp_func = interp1d(summed_spectra['wav'], summed_spectra['flux'])
flux_sum = interp_func(const_grid)

const_gridd_dnr = np.arange(np.min(donor_spectra['wav']), np.max(donor_spectra['wav']), grid_size)
interp_func_dnr = interp1d(donor_spectra['wav'], donor_spectra['flux'])
flux_dnr = interp_func_dnr(const_gridd_dnr)

## Instrumental broadening ##
from astropy.convolution import Gaussian1DKernel, convolve

FWHM = 6 #AA approx. 

sd = FWHM / (grid_size * 2*np.sqrt(2 * np.log(2))) # want sigma instead of FMHM
g = Gaussian1DKernel(stddev=sd)
z_sum = convolve(flux_sum, g)
z_donor = convolve(flux_dnr, g)

## Rotational broadening ##
from PyAstronomy import pyasl
P = 0.11414278 #days 
eps = 0.5 # linear limb-darkening coeff
rv = (2 * np.pi * 0.23 * 695700)/(P*86400)  # km/s, projected rotational velocity

rflux_sum = pyasl.rotBroad(const_grid, z_sum, epsilon = eps, vsini = rv)
rflux_dnr = pyasl.rotBroad(const_gridd_dnr, z_donor, epsilon = eps, vsini = rv)

## Normalization ##
rolling_sum = pd.Series(rflux_sum).rolling(int(300/grid_size), min_periods = 1, center = True).median()
rolling_dnr = pd.Series(rflux_dnr).rolling(int(300/grid_size), min_periods = 1, center = True).median()

spec_WD_donor = pd.DataFrame(data = {'wavelength': const_grid, 'flux': rflux_sum/rolling_sum})
spec_donor = pd.DataFrame(data = {'wavelength': const_gridd_dnr, 'flux': rflux_dnr/rolling_dnr})

spec_WD_donor.to_csv(save_loc + 'spec_WD_donor.csv', index = False)
spec_donor.to_csv(save_loc + 'spec_donor.csv', index = False)

#%% Making plots 

## Full spectra (broadened but not normalized) ##

plt.plot(WD['wav'], WD['flux'], lw = 0.3, label = 'WD, T = 20000K, logg = 8.0')
plt.plot(const_gridd_dnr, rflux_dnr, lw = 0.3,  label = 'donor, T = 5000K, logg = 4.5')
plt.plot(const_grid, rflux_sum,  lw = 0.3, label = 'WD + donor')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux (erg s$^{-1}$cm$^{-2}$$\AA$$^{-1}$)')
plt.legend()
plt.savefig(save_loc + 'full_spectra_brdnd.png', dpi = 300, bbox_inches = 'tight')
plt.close()

## Comparison ##

plt.plot(spec_WD_donor[(spec_WD_donor['wavelength'] > 1000) & (spec_WD_donor['wavelength'] < 15000)]['wavelength'], 
         spec_WD_donor[(spec_WD_donor['wavelength'] > 1000) & (spec_WD_donor['wavelength'] < 15000)]['flux'], 
         lw = 0.3, label = 'WD + donor, brdnd + norm')
plt.plot(spec_donor[(spec_donor['wavelength'] > 1000) & (spec_donor['wavelength'] < 15000)]['wavelength'], 
         spec_donor[(spec_donor['wavelength'] > 1000) & (spec_donor['wavelength'] < 15000)]['flux'], 
         lw = 0.3, label = 'donor, brdnd + norm')
plt.legend()
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.savefig(save_loc + 'both_brdnd_norm_full.png', dpi = 300, bbox_inches = 'tight')
plt.close()

plt.plot(spec_WD_donor[(spec_WD_donor['wavelength'] > 3500) & (spec_WD_donor['wavelength'] < 10000)]['wavelength'], 
         spec_WD_donor[(spec_WD_donor['wavelength'] > 3500) & (spec_WD_donor['wavelength'] < 10000)]['flux'], 
         lw = 0.3, label = 'WD + donor, brdnd + norm')
plt.plot(spec_donor[(spec_donor['wavelength'] > 3500) & (spec_donor['wavelength'] < 10000)]['wavelength'], 
         spec_donor[(spec_donor['wavelength'] > 3500) & (spec_donor['wavelength'] < 10000)]['flux'], 
         lw = 0.3, label = 'donor, brdnd + norm')
plt.legend()
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.savefig(save_loc + 'both_brdnd_norm.png', dpi = 300, bbox_inches = 'tight')
plt.close()

# Ca line

plt.plot(spec_WD_donor[(spec_WD_donor['wavelength'] > 3900) & (spec_WD_donor['wavelength'] < 4000)]['wavelength'], 
         spec_WD_donor[(spec_WD_donor['wavelength'] > 3900) & (spec_WD_donor['wavelength'] < 4000)]['flux'], 
         lw = 0.3, label = 'WD + donor, brdnd + norm')
plt.plot(spec_donor[(spec_donor['wavelength'] > 3900) & (spec_donor['wavelength'] < 4000)]['wavelength'], 
         spec_donor[(spec_donor['wavelength'] > 3900) & (spec_donor['wavelength'] < 4000)]['flux'], 
         lw = 0.3, label = 'donor, brdnd + norm')
plt.legend()
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.savefig(save_loc + 'both_brdnd_norm_Ca.png', dpi = 300, bbox_inches = 'tight')
plt.close()

## Individual  ##

plt.plot(spec_WD_donor[(spec_WD_donor['wavelength'] > 3500) & (spec_WD_donor['wavelength'] < 10000)]['wavelength'], 
         spec_WD_donor[(spec_WD_donor['wavelength'] > 3500) & (spec_WD_donor['wavelength'] < 10000)]['flux'], 
         lw = 0.3, label = 'WD + donor, brd + norm')
plt.legend()
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.savefig(save_loc + 'sum_brdnd_norm.png', dpi = 300, bbox_inches = 'tight')
plt.close()

plt.plot(spec_donor[(spec_donor['wavelength'] > 3500) & (spec_donor['wavelength'] < 10000)]['wavelength'], 
         spec_donor[(spec_donor['wavelength'] > 3500) & (spec_donor['wavelength'] < 10000)]['flux'], 
         lw = 0.3, label = 'donor, brdnd + norm')
plt.legend()
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.savefig(save_loc + 'donor_brdnd_norm.png', dpi = 300, bbox_inches = 'tight')

