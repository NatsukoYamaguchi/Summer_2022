#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 09:10:29 2022

@author: natsukoyamaguchi

Finds best fit parameters of lcurve to model the observed light curve using the maximum likelihood method.

"""
from lcurve_functions import lcurve
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize 
from pathlib import Path

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/lcurve_fitting'
Path(folder_dir + '/' + 'ml_fitting').mkdir(parents=True, exist_ok=True)

# file_name = 'IC_2_zr'
file_name = 'tm_t1_25000'

observed = pd.read_csv(folder_dir +'/observed/new_data.csv')
# Depending on the filter (r/g), make sure to set ldc, gc, and wavelength to the appropriate values. 

x = observed['phase']
y = observed['flux']
yerr = observed['flux_err']

# Maximum likelihood function (chi square)

def log_likelihood(theta, x, y, yerr):
    iangle, t1, rdisc2, temp_disc, height_disc, texp_disc, t0 = theta 
    model = lcurve(iangle = iangle, t1 = t1, rdisc2 = rdisc2, temp_disc = temp_disc, height_disc = height_disc, texp_disc = texp_disc, t0 = t0)
    func = interp1d(model['phase'], model['flux_norm'], fill_value='extrapolate') 
    model_y = func(x) # interpolate to the observed phases 
    sigma2 = yerr**2 
    return -0.5 * np.sum((y - model_y)**2 / sigma2 + np.log(sigma2)) # ln of eqn (9) in Hogg et al. 

bnds = ((60, 90), (25000, 35000), (0.2, 0.5), (2000, 7000), (0.05, 0.5), (-3, -0.5), (56489.31, 56489.33)) # bounds on each parameter
initial = np.array([80, 28000, 0.3, 3000, 0.1, -1.5, 56489.32])

nll = lambda *args: -log_likelihood(*args) 

# minimize the negative of the log likelihood func (maximize the log likelihood)
soln = minimize(nll, initial, args=(x, y, yerr), method = 'Nelder-Mead', bounds = bnds)

iangle_ml, t1_ml, rdisc2_ml, temp_disc_ml, height_disc_ml, texp_disc_ml, t0_ml = soln.x
dat_ml = lcurve(iangle = iangle_ml, t1 = t1_ml, rdisc2 = rdisc2_ml, temp_disc = temp_disc_ml, height_disc = height_disc_ml, texp_disc = texp_disc_ml, t0 = t0_ml)

#%% 

# Plot the results 

plt.errorbar(x, y, yerr, fmt=".k", capsize=0, label = "New data")
plt.plot(dat_ml['phase'], dat_ml['flux_norm'] , label="ML")
plt.legend()
plt.xlabel("phase")
plt.ylabel("flux")

plt.text(1.1, 1, 'ML parameters: \n iangle = {0:.3f}'.format(iangle_ml) + ' \n t1 = {0:.3f}'.format(t1_ml) + ' \n rdisc2 = {0:.3f}'.format(rdisc2_ml) + ' \n temp_disc = {0:.3f}'.format(temp_disc_ml) + ' \n height_disc = {0:.3f}'.format(height_disc_ml) + ' \n texp_disc = {0:.3f}'.format(texp_disc_ml) +  ' \n t0 = {0:.3f}'.format(t0_ml), fontsize = 10)
plt.text(1.1, 0.7, 'Initial values: \n iangle = {0:.3f}'.format(initial[0]) + ' \n t1 = {0:.3f}'.format(initial[1]) + ' \n rdisc2 = {0:.3f}'.format(initial[2]) + ' \n temp_disc = {0:.3f}'.format(initial[3]) + ' \n height_disc = {0:.3f}'.format(initial[4])  + ' \n texp_disc = {0:.3f}'.format(initial[5]) +  ' \n t0 = {0:.3f}'.format(initial[6]), fontsize = 10)

plt.savefig(folder_dir + '/' +  'ml_fitting' + '/' + file_name + '.png', dpi = 300, bbox_inches =  "tight")

# Save ML parameters

params = pd.Series(data = {'iangle': iangle_ml, 't1': t1_ml, 'rdisc2': rdisc2_ml, 'temp_disc': temp_disc_ml, 'height_disc': height_disc_ml, 'texp_disc':texp_disc_ml, 't0':t0_ml})
params.to_csv(folder_dir + '/' +  'ml_fitting' + '/' + '/params_' + file_name + '.csv')
