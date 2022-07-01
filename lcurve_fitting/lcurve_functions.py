#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:40:15 2022

@author: natsukoyamaguchi
"""

import pexpect
import numpy as np
import pandas as pd
from pathlib import Path

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/lcurve_fitting'
Path(folder_dir + '/' + 'lroche_files').mkdir(parents=True, exist_ok=True)

param_file_dir = folder_dir + '/' + 'lroche_files' + '/param_file'
dat_file_dir = folder_dir + '/' + 'lroche_files' + '/dat'

def create_param_file(
        save_loc = param_file_dir,
        q=0.2096774194,
        iangle=79.76,
        r1=0.005,
        r2=-1,
        cphi3=0.015,
        cphi4=0.017,
        t1=10000,
        t2=4726,
        spin1=1,
        spin2=1,
        ldc1_1=0.81760109999999997,
        ldc1_2=-0.66854585,
        ldc1_3=0.61190517,
        ldc1_4=-0.23623394,
        ldc2_1=0.6272,   # r: 0.6399, g: 0.6272
        ldc2_2=-0.8987,  # r: -0.6048, g: -0.8987	
        ldc2_3=1.8536,   # r: 1.3699, g: 1.8536
        ldc2_4=-0.6704,  # r: -0.5698, g:-0.6704
        velocity_scale=323,
        beam_factor1=1.31312039536,
        beam_factor2=1.31312039536,
        deltat=0,
        t0=56489.38307139,  #56489.326 + P/2 (Birth of the ELMs paper has phases off by 0.5)
        period=0.11414278,
        gravity_dark1=0.5,
        gravity_dark2=0.7501, # r: 0.5424, g: 0.7501
        absorb=1,
        slope=0,
        cube=0,
        quad=0,
        rdisc1=-1,
        rdisc2=0.3,
        height_disc=0.2,
        beta_disc=1.5,
        temp_disc=3200,
        texp_disc=-1.8,
        temp_edge=4000,
        absorb_edge=1,
        lin_limb_disc=0.3,
        quad_limb_disc=0,
        radius_spot=0.34999999999999998,
        length_spot=0.02,
        height_spot=0.050000000000000003,
        expon_spot=0,
        epow_spot=1,
        angle_spot=140,
        yaw_spot=0,
        temp_spot=15000,
        tilt_spot=90,
        cfrac_spot=0.20000000000000001,
        delta_phase=3e-08,
        nlat1f=50,
        nlat2f=210,
        nlat1c=50,
        nlat2c=170,
        npole=1,
        nlatfill=3,
        nlngfill=2,
        lfudge=0.07,
        llo=0,
        lhi=-50,
        phase1=0.018,
        phase2=0.482,
        wavelength=600,
        roche1=0,
        roche2=1,
        eclipse1=1,
        eclipse2=1,
        glens1=1,
        tperiod=0.1141427,
        gdark_bolom1=0,
        gdark_bolom2=0,
        mucrit1=0,
        mucrit2=0,
        limb1='Claret',
        limb2='Claret',
        use_radii=1,
        mirror=0,
        add_disc=1,
        nrad=40,
        opaque=1,
        add_spot=0,
        nspot=100,
        iscale=0): 
    
    '''
    Creates a parameter file to be used as the model for lroche. 

    Parameters
    ----------
    save_loc : string
        Directory where the output parameter file should be saved. If a file with the same name already exists, it will be overwritten. Otherwise, a new file will be created.
    **kwargs: parameters for lroche. All floats except for limb1 and limb2 which are strings. For more details, see documentation for lroche.
    
    Returns
    -------
    None. Creates a parameter file in specified directory.

    '''

    
    with open(save_loc, 'w') as f:
        f.write('q                    = ' + str(q) + ' 0.02 0.01 1\n')
        f.write('iangle               = ' + str(iangle) + ' 0.2 0.03 1\n')
        f.write('r1                   = ' + str(r1) + ' 0.0002 0.0001 1\n')
        f.write('r2                   = ' + str(r2) + ' 0.005 0.0002 1\n')
        f.write('cphi3                = ' + str(cphi3) + ' 0.001 0.001 0\n')
        f.write('cphi4                = ' + str(cphi4) + ' 0.001 0.001 0\n')
        f.write('t1                   = ' + str(t1) + ' 50 10 1\n')
        f.write('t2                   = ' + str(t2) + ' 10 10 1\n')
        f.write('spin1                = ' + str(spin1) + ' 0.001 0.001 0\n')
        f.write('spin2                = ' + str(spin2) + ' 0.001 0.001 0\n')
        f.write('ldc1_1               = ' + str(ldc1_1) + ' 0.01 0.01 0\n')
        f.write('ldc1_2               = ' + str(ldc1_2) + ' 0.01 0.01 0\n')
        f.write('ldc1_3               = ' + str(ldc1_3) + ' 0.01 0.01 0\n')
        f.write('ldc1_4               = ' + str(ldc1_4) + ' 0.01 0.01 0\n')
        f.write('ldc2_1               = ' + str(ldc2_1) + ' 0.01 0.01 0\n')
        f.write('ldc2_2               = ' + str(ldc2_2) + ' 0.01 0.01 0\n')
        f.write('ldc2_3               = ' + str(ldc2_3) + ' 0.01 0.01 0\n')
        f.write('ldc2_4               = ' + str(ldc2_4) + ' 0.01 0.01 0\n')
        f.write('velocity_scale       = ' + str(velocity_scale) + ' 1 1 1\n')
        f.write('beam_factor1         = ' + str(beam_factor1) + ' 0.1 0.02 0\n')
        f.write('beam_factor2         = ' + str(beam_factor2) + ' 0.1 0.002 1\n')
        f.write('deltat               = ' + str(deltat) + ' 0.001 0.001 0\n')
        f.write('t0                   = ' + str(t0) + ' 0.0001 1e-05 1\n')
        f.write('period               = ' + str(period) + ' 1e-06 1e-06 0\n')
        f.write('gravity_dark1        = ' + str(gravity_dark1) + ' 0.0001 0.0001 0\n')
        f.write('gravity_dark2        = ' + str(gravity_dark2) + ' 0.0001 0.0001 0\n')
        f.write('absorb               = ' + str(absorb) + ' 0.001 0.001 0\n')
        f.write('slope                = ' + str(slope) + ' 0.01 1e-05 0\n')
        f.write('cube                 = ' + str(cube) + ' 0.01 1e-05 0\n')
        f.write('quad                 = ' + str(quad) + ' 0.01 1e-05 0\n')
        f.write('rdisc1               = ' + str(rdisc1) + ' 0.01 0.0001 0\n')
        f.write('rdisc2               = ' + str(rdisc2) + ' 0.01 0.02 0\n')
        f.write('height_disc          = ' + str(height_disc) + ' 0.01 1e-05 0\n')
        f.write('beta_disc            = ' + str(beta_disc) + ' 0.01 1e-05 0\n')
        f.write('temp_disc            = ' + str(temp_disc) + ' 50 40 0\n')
        f.write('texp_disc            = ' + str(texp_disc) + ' 0.2 0.001 0\n')
        f.write('temp_edge            = ' + str(temp_edge) + ' 0.001 0.001 0\n')
        f.write('absorb_edge          = ' + str(absorb_edge) + ' 0.001 0.001 0\n')
        f.write('lin_limb_disc        = ' + str(lin_limb_disc) + ' 0.02 0.0001 0\n')
        f.write('quad_limb_disc       = ' + str(quad_limb_disc) + ' 0.02 0.0001 0\n')
        f.write('radius_spot          = ' + str(radius_spot) + ' 0.005 0.01 0\n')
        f.write('length_spot          = ' + str(length_spot) + ' 0.002 0.005 0\n')
        f.write('height_spot          = ' + str(height_spot) + ' 0.005 1e-05 0\n')
        f.write('expon_spot           = ' + str(expon_spot) + ' 0.2 0.1 0\n')
        f.write('epow_spot            = ' + str(epow_spot) + ' 0.1 0.1 0\n')
        f.write('angle_spot           = ' + str(angle_spot) + ' 2 2 0\n')
        f.write('yaw_spot             = ' + str(yaw_spot) + ' 2 2 0\n')
        f.write('temp_spot            = ' + str(temp_spot) + ' 400 200 0\n')
        f.write('tilt_spot            = ' + str(tilt_spot) + ' 5 2 0\n')
        f.write('cfrac_spot           = ' + str(cfrac_spot) + ' 0.02 0.008 0\n')
        f.write('delta_phase          = ' + str(delta_phase) + '\n')
        f.write('nlat1f               = ' + str(nlat1f) + '\n')
        f.write('nlat2f               = ' + str(nlat2f) + '\n')
        f.write('nlat1c               = ' + str(nlat1c) + '\n')
        f.write('nlat2c               = ' + str(nlat2c) + '\n')
        f.write('npole                = ' + str(npole) + '\n')
        f.write('nlatfill             = ' + str(nlatfill) + '\n')
        f.write('nlngfill             = ' + str(nlngfill) + '\n')
        f.write('lfudge               = ' + str(lfudge) + '\n')
        f.write('llo                  = ' + str(llo) + '\n')
        f.write('lhi                  = ' + str(lhi) + '\n')
        f.write('phase1               = ' + str(phase1) + '\n')
        f.write('phase2               = ' + str(phase2) + '\n')
        f.write('wavelength           = ' + str(wavelength) + '\n')
        f.write('roche1               = ' + str(roche1) + '\n')
        f.write('roche2               = ' + str(roche2) + '\n')
        f.write('eclipse1             = ' + str(eclipse1) + '\n')
        f.write('eclipse2             = ' + str(eclipse2) + '\n')
        f.write('glens1               = ' + str(glens1) + '\n')
        f.write('tperiod              = ' + str(tperiod) + '\n')
        f.write('gdark_bolom1         = ' + str(gdark_bolom1) + '\n')
        f.write('gdark_bolom2         = ' + str(gdark_bolom2) + '\n')
        f.write('mucrit1              = ' + str(mucrit1) + '\n')
        f.write('mucrit2              = ' + str(mucrit2) + '\n')
        f.write('limb1                = ' + limb1 + '\n')
        f.write('limb2                = ' + limb2 + '\n')
        f.write('use_radii            = ' + str(use_radii) + '\n')
        f.write('mirror               = ' + str(mirror) + '\n')
        f.write('add_disc             = ' + str(add_disc) + '\n')
        f.write('nrad                 = ' + str(nrad) + '\n')
        f.write('opaque               = ' + str(opaque) + '\n')
        f.write('add_spot             = ' + str(add_spot) + '\n')
        f.write('nspot                = ' + str(nspot) + '\n')
        f.write('iscale               = ' + str(iscale) + '')		
        
def run_lroche(time1 = 59619.75221737, 
               time2 = 59619.87119743, 
               ntime = 5000,
               expose = 0,
               param_file_loc = param_file_dir, save_loc = dat_file_dir):
    
    '''
    Runs lroche in trm_software. 

    Parameters
    ----------
    param_file_loc : string
        Directory of the parameter file. 
    save_loc : string
        Directory of where the output data file will be saved. 

    Returns
    -------
    None. Just saves the output as a text file. 

    '''
    
    p = pexpect.spawn('/Users/natsukoyamaguchi/trm_software/bin/lcurve/lroche')
    p.expect(':') # MODEL - model file of parameter values
    p.sendline(param_file_loc)
    p.expect(':') # DATA - data file
    p.sendline('none')
    p.expect(':') # TIME1 - first time to compute
    p.sendline(str(time1))
    p.expect(':') # TIME2 - last time to compute  
    p.sendline(str(time2))
    p.expect(':') # NTIME - number of times to compute 
    p.sendline(str(ntime))
    p.expect(':') # EXPOSE - exposure time (same units as ephemeris)
    p.sendline(str(expose))
    p.expect(':') # NDIVIDE - ndivide factor to use 
    p.sendline('1')
    p.expect(':') # NOISE - RMS noise to add (after scaling) 
    p.sendline('0')
    p.expect(':') # SEED - random number seed
    p.sendline('57565')
    p.expect(':') # NFILE - number of files to generate
    p.sendline('1')
    p.expect(':') # OUTPUT - file/root to save data
    p.sendline(save_loc)
    p.expect(':') # DEVICE - plot device
    p.sendline('none') # '/xs' to display the plot in XQuartz 
    p.expect(':') # SSFAC - global scale factor
    p.sendline('2504079.7566472059')
    p.wait()
    

def read_lightcurve(save_loc = dat_file_dir, period = 0.11414278):
    '''
    Reads in lightcurve generated by lroche

    Parameters
    ----------
    save_loc : string
        Location of the dat file outputted by lroche. 

    Returns
    -------
    df : pandas dataframe
        The columns are midexp_times, phase, flux, flux_norm (sorted by phase)

    '''
    dat = np.loadtxt(save_loc)
    
    
    phase = dat[:,0]  % period / period
    
    # df = pd.DataFrame({'midexp_times': dat[:,0], 'exp_times':dat[:,1], 'phase': phase, 'flux': dat[:,2], 'flux_norm': dat[:,2]/np.median(dat[:,2]), 'flux_uncert':dat[:,3], 
    #                'weight': dat[:,4], 'int': dat[:,5]})
    
    df = pd.DataFrame({'midexp_times': dat[:,0], 'phase': phase, 'flux': dat[:,2], 'flux_norm': dat[:,2]/np.median(dat[:,2])})
    
    df2 = df.sort_values('phase')
    
    return df2

def lcurve(**params):

    '''
    
    Parameters
    ----------
    params : kwargs for parameters for lroche 
        Default values will be set for any that are unspecified.
        For a list of possible parameters and their default values, look at docstring for create_param_file(). 

    Returns
    -------
    lcurve_df : pandas dataframe 
        Data frame with the output from lroche given the inputted parameters 
        The columns are midexp_times, phase, flux, flux_norm (sorted by phase)
    '''
    
    create_param_file(**params)
    
    run_lroche()
 
    data = read_lightcurve()
    
    return data