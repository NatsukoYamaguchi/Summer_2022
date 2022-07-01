#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:46:05 2022

@author: natsukoyamaguchi

Get RV vs phase (+ sinusoidal fit) from the spectra (using CCF)

"""

# base_settings -> normalization -> skyline -> broadening -> RV

from base_settings import *
from normalization import *
from skyline import *
from temp_broadening import *

broadening = 'yes' # should be left as 'yes' ('no' means not accounting for instantaneous or rotational broadening)

#%% Get the normalized spectrum 

save_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/outputs'

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
            dat_i = data = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + files[i])
        elif norm == 'rolling':
            dat_i = data = pd.read_csv(base_directory + '/spectra_norm/' + norm + '/' + str(window) + '/'  + files[i])
        dat.append(dat_i)

#%% Download and normalize template spectra

if col == 'r':
    x_lolim = 5400 # wavelength range used to match with template 
    x_uplim = 6600
elif col == 'b':
    x_lolim = 3800
    x_uplim = 5200
    
from astropy.io import fits

if broadening == 'no': 
    # BOSZ Kurucz Model - T_eff = 5000K, log g = 4.5, 0, 0, 0, Inst broad = 2000
    hdul = fits.open(base_directory + 'templates/amp00cp00op00t5000g45v20modrt0b2000rs.fits')
    template = hdul[1].data

    # .byteswap().newbyteorder()
    # for information, go to https://github.com/astropy/astropy/issues/1156
    x_temp = template.Wavelength.byteswap().newbyteorder()
    y_temp = template.SpecificIntensity.byteswap().newbyteorder()
    
    # Isolating relevant wavelength range and fit a polynomial
    
    x_temp_slice = x_temp[np.logical_and(x_temp >= x_lolim, x_temp <= x_uplim)]
    y_temp_slice = y_temp[np.logical_and(x_temp >= x_lolim, x_temp <= x_uplim)]
    
    fit_poly = np.polyfit(x_temp_slice, y_temp_slice, 2)
    polyfit_norm = fit_poly[0]*x_temp_slice**2 + fit_poly[1]*x_temp_slice + fit_poly[2]
    
    flux_norm_temp = y_temp_slice/polyfit_norm
     
    # *log* of wavelength (the dopper shift is constant on log scale for the CCF)
    temp_spectrum = pd.DataFrame(data = {'wavelength': np.log(x_temp_slice), 'flux': flux_norm_temp})
    tag = 'no_inst_rot'

elif broadening == 'yes': 
    if norm == 'poly':
        temp_spectrum = pd.read_csv(save_dir + '/temp_inst_rot_' + norm + '_' + col + '.csv')
        tag = 'inst_rot' 
    elif norm == 'rolling':
        temp_spectrum = pd.read_csv(save_dir + '/temp_inst_rot_' + norm + '_' + str(window) + '_' + col + '.csv')
        tag = 'inst_rot' 

#%% Cross-correlation 

from scipy.interpolate import interp1d

v = np.arange(0.001, 600, 10) # km/s (not starting at zero to avoid zero error)
shift = np.sort(np.concatenate((-np.log(1 + v/(3e5)), np.log(1 + v/(3e5)))))

def cov(X, Y):
    ''' Covariance of X and Y (must be of the same size) '''
    X_bar = np.mean(X)
    Y_bar = np.mean(Y)
    return np.mean((X - X_bar)*(Y - Y_bar))

def CCF(X, Y):
    ''' Cross-correlation function as defined in Appendix A of Zhang et al. 2021 '''
    return cov(X,Y)/(np.sqrt(np.var(X)*np.var(Y)))

shift_arr = np.zeros(len(dat))

for n in range(len(dat)):
    # x_dat_n = dat[n][dat[n]['flux_norm'] >= 0]['wavelength'] # Mask away points with negative fluxes 
    # y_dat_n = dat[n][dat[n]['flux_norm'] >= 0]['flux_norm']
    x_dat_n = dat[n]['wavelength']  
    y_dat_n = dat[n]['flux_norm']
    x_dat = x_dat_n[(x_dat_n >= x_lolim + 200) & (x_dat_n <= x_uplim - 200)]
    y_dat = y_dat_n[(x_dat_n >= x_lolim + 200) & (x_dat_n <= x_uplim - 200)]
    
    # *log* of wavelength (the dopper shift is constant on log scale for the CCF)
    dat_spectrum = pd.DataFrame(data = {'wavelength': np.log(x_dat), 'flux': y_dat})
    
    CCF_arr = np.zeros(len(shift))
    
    for i in range(len(shift)):
    
        temp_shifted = pd.DataFrame(data = {'wavelength': temp_spectrum['wavelength'] + shift[i], 'flux': temp_spectrum['flux']})
        interp_func = interp1d(temp_shifted['wavelength'], temp_shifted['flux'])
        
        temp_flux = interp_func(dat_spectrum['wavelength'])
        
        CCF_arr[i] = CCF(dat_spectrum['flux'], temp_flux)
    
    def func(x):
        
        interp_func = interp1d(temp_spectrum['wavelength'] + x, temp_spectrum['flux'])
        return 1/CCF(dat_spectrum['flux'], interp_func(dat_spectrum['wavelength']))
    
    # use the shift that maximized CCF from the range we tried as an initial guess to find the precise shift 
    res = optimize.minimize(func, x0 = shift[np.argmax(CCF_arr)], method='nelder-mead')
    shift_arr[n] = res.x
   
#%% Plot RV vs orbital phase

RV = (np.exp(shift_arr) -1) * 3e5

x = phase_list[phase_list.argsort()]
y = RV[phase_list.argsort()]

def sin_fit(z, a, b, c):
    return a * np.sin(b * z) + c

params, params_covariance = optimize.curve_fit(sin_fit, x, y,
                                               p0=[350, 6.28, -100])

# plt.xlabel('Orbital Phase')
# plt.ylabel('RV (km/s)')
# plt.scatter(x, y, s = 3)
# plt.plot(x, sin_fit(x, params[0], params[1], params[2]))
# plt.scatter(-x, y, s = 3)
# plt.plot(-x, cos_fit(x, params[0], params[1], params[2]), label='Fitted function')

#%% Barycenter correction (BC)

from astropy.time import Time 
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

# modify ra and dec format 
ra_ha = []
dec_deg = []

for i in range(len(ra)):
    ra_ha.append((float(ra[i].split(':')[0]) + float(ra[i].split(':')[1])/60 + float(ra[i].split(':')[2])/3600) * u.hourangle)
     
    if dec[i][0] == '+':
        dec_ = float(dec[i].split(':')[0].replace('+', ''))
    else:
        dec_ = float(dec[i].split(':')[0])

    dec_deg.append((dec_ * u.deg + float(dec[i].split(':')[1]) * u.arcmin + float(dec[i].split(':')[2]) * u.arcsec).to(u.deg))

# Calculate barycentre correction at each observation 
barycorr = np.zeros(len(dat))
for n in range(len(dat)):
    time = Time(jdmid[n], format = 'jd')
    # keck = EarthLocation.of_site('Keck') # needs internet
    keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)  
    sc = SkyCoord(ra= ra_ha[n], dec= dec_deg[n])
    barycorr[n] = sc.radial_velocity_correction(obstime=time, location=keck).to(u.km/u.s).value

print('average barycentre correction:', np.mean(barycorr))

params_BC, params_covariance_BC = optimize.curve_fit(sin_fit, phase_list[phase_list.argsort()], RV[phase_list.argsort()] + barycorr,
                                               p0=[350, 6.28, -100])

# plt.xlabel('Orbital Phase')
# plt.ylabel('RV (km/s)')
# plt.scatter(phase_list[phase_list.argsort()], RV[phase_list.argsort()], s = 3)
# plt.plot(phase_list[phase_list.argsort()], sin_fit(phase_list[phase_list.argsort()], params[0], params[1], params[2]),
#           label = 'without BC correction', lw = 0.8)
# plt.scatter(phase_list[phase_list.argsort()], RV[phase_list.argsort()] + barycorr, s = 3)
# plt.plot(phase_list[phase_list.argsort()], sin_fit(phase_list[phase_list.argsort()], params_BC[0], params_BC[1], params_BC[2]),
#           label='with BC correction', lw = 0.8)
# plt.legend()

#%%

df = pd.DataFrame(data = {'phase': phase_list[phase_list.argsort()], 'RV': RV[phase_list.argsort()] + barycorr, 'fit': sin_fit(phase_list[phase_list.argsort()], params_BC[0], params_BC[1], params_BC[2])})
if norm == 'poly':
    fn = 'RV_' + 'BC_' + tag + '_' + norm + '_' + col + '.csv'
elif norm == 'rolling':
    fn = 'RV_' + 'BC_' + tag + '_' + norm + '_' + str(window) + '_' + col + '.csv'
df.to_csv(base_directory + '/outputs/' + fn, index = False)