#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:15:41 2022

@author: natsukoyamaguchi

Get time information from the comments of the lris files 
+ normalize the spectra (either by polynomial or rolling median fit the the continuum)

"""
from base_settings import *

Path(base_directory + '/spectra_norm/poly').mkdir(parents=True, exist_ok=True)
Path(base_directory + '/spectra_norm/rolling').mkdir(parents=True, exist_ok=True)
Path(base_directory + '/spectra_norm/rolling' + '/' + str(window)).mkdir(parents=True, exist_ok=True)

#%% Get the lris spectrum 

files = sorted(listdir(base_directory + '/lris_p274')) # alphabetical order

if col == 'r':
    label = 'celsr'
    x_lolim = 5400 # wavelength range used to match with template 
    x_uplim = 6600
elif col == 'b':
    label = 'celsb'
    x_lolim = 3800
    x_uplim = 5200

dat_nan = []
for i in range(len(files)):
    if files[i][0:5] == label:
        dat_nan.append(np.loadtxt(base_directory + '/lris_p274' + '/' + files[i])) 

dat = []

for i in range(len(dat_nan)):
    dat_i = dat_nan[i]
    dat.append(dat_i[~np.isnan(dat_i).any(axis=1)]) # erase all rows with NaNs
    
# Get the JDMID from comments
    
jdmid = []
ra = []
dec = []

for i in range(len(files)):
    with open(base_directory + '/lris_p274' + '/' + files[i]) as f:
        lines = f.readlines()[:122]
        JDMID = lines[22]
        if label == 'celsr':
            RA = lines[30]
            DEC = lines[31]
        elif label == 'celsb':
            RA = lines[32]
            DEC = lines[33]
        
    if files[i][0:5] == label:
        jdmid.append(float(JDMID.replace(" ", "")[7:20])) #stores just the JDMID as float
        ra.append(RA.replace(" ", "")[5:16])
        dec.append(DEC.replace(" ", "")[6:17])

phase_list = np.zeros(len(jdmid))

for i in range(len(jdmid)):
    phase = jdmid[i] % (P) / (P) # convert dates into phase (2.74 hours)
    phase_list[i] = phase

# Save the times
    
df = pd.DataFrame(data={'jdmid': jdmid, 'phase': phase_list, 'ra': ra, 'dec': dec})
df.to_csv(base_directory + '/spectra_norm' + '/' + col + '_' + 'times' + '.csv', index = False)

# Normalize spectra by fitting a polynomial or rolling median to the continuum
             
if norm == 'poly':
    
    for i in range(len(dat)):
        x = dat[i][:,0]  # wavelength
        y = dat[i][:,1]  # flux 
        y_2 = dat[i][:,2]  # sky flux 
        z = np.polyfit(x, y, 2) # spectrum continuum
        z_2 = np.polyfit(x, y_2, 2) # sky continuum
        y_fit = z[0] * x**2 + z[1] * x + z[2]
        y_fit_2 = z_2[0] * x**2 + z_2[1] * x + z_2[2]
        df = pd.DataFrame(data={'wavelength': x, 'flux': y, 'sky_flux': y_2, 'flux_norm': y/y_fit, 'sky_flux_norm': y_2/y_fit_2})
        df.to_csv(base_directory + '/spectra_norm/poly' + '/' + col + '_' + str(jdmid[i]) + '.csv', index = False)
    
elif norm == 'rolling':

    for i in range(len(dat)):
        x = dat[i][:,0] # wavelength   
        y = dat[i][:,1] # flux 
        y_2 = dat[i][:,2] #sky_flux
        fine_grid = np.arange(np.min(x), np.max(x), 0.14) # interpolate to finer grid ~ same as template 
        interp_func = interp1d(x, y) 
        interp_func_2 = interp1d(x, y_2)
        rolling = pd.Series(interp_func(fine_grid)).rolling(window, min_periods = 1, center = True).median()
        rolling_2 = pd.Series(interp_func_2(fine_grid)).rolling(window, min_periods = 1, center = True).median()
        interp_func_back = interp1d(fine_grid, rolling, fill_value="extrapolate")
        interp_func_back_2 = interp1d(fine_grid, rolling_2, fill_value="extrapolate")
        cont = y/interp_func_back(x)
        cont_2 = y_2/interp_func_back_2(x) # back to original grid 
        df = pd.DataFrame(data={'wavelength': x, 'flux': y, 'sky_flux': y_2, 'flux_norm': cont, 'sky_flux_norm': cont_2})
        df.to_csv(base_directory + '/spectra_norm/rolling' + '/' + str(window) + '/' + col + '_' +  str(jdmid[i]) + '.csv', index = False)

#%%

# plt.close()

# i = 2
# x = dat[i][:,0]  
# # y = dat[i][:,1]
# y = dat[i][:,2]  
# df_temp = pd.DataFrame(data = {'wavelength': x, 'flux': y, 'sky_flux': y})
# rolling = df_temp['flux']/(df_temp['flux'].rolling(100, min_periods = 1).median())
# rolling2 = df_temp['sky_flux']/(df_temp['sky_flux'].rolling(100, min_periods = 1).median())    


# plt.plot(x, y, lw = 0.3)
# plt.plot(x, df_temp['flux'].rolling(100, min_periods = 1).median())


# diff = np.zeros(len(df_temp['wavelength'])-1)
# for j in range(len(df_temp['wavelength'])-1):
#     diff[j] = df_temp['wavelength'][j+1] - df_temp['wavelength'][j]


# wavelength_grid = np.arange(np.min(x), np.max(x), 0.14)
# interp_func = interp1d(x, y)

# #%%
# plt.plot(wavelength_grid, interp_func(wavelength_grid), lw = 0.3)
# series = pd.Series(interp_func(wavelength_grid))
# plt.plot(wavelength_grid, series.rolling(1000, min_periods = 1).median())

# plt.show()

# #%% 

# interp_func_2 = interp1d(wavelength_grid, series.rolling(500, min_periods = 1).median(), fill_value="extrapolate")
# cont = interp_func_2(x)

# plt.plot(x, y, lw = 0.3)
# plt.plot(x, cont)

# #%%

# plt.plot(x, y/cont)


