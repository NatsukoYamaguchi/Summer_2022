#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 15:33:50 2022

@author: natsukoyamaguchi
"""

import numpy as np 

def get_koester_spectrum_cgs_units(path= 'data/koester_wd_models/koester_wd_20kK_logg8.dat', d_pc = 1000, 
    mass_assume = 0.6):
    '''
    Spectrum depends on Teff and R, but Koester models are in terms of Teff and logg. 
    This calculates the appropriate R and scales the spectrum accordingly.
    output in erg/s/cm^2/Angstrom
    '''
    from astropy import units as u, constants as c
    
    dat = np.genfromtxt(path)
    lines = open(path).readlines()[2]
    logg = float(lines.split('logg = ')[1].split('log')[0]) 
    
    # Mass-radius relation from Bedard et al. 
    m_grid = np.array([0.2, 0.25, 0.30, 0.35, 0.40 , 0.45, 0.50, 0.55, 0.60, 0.65,0.70, 0.75, 0.80, 0.85, 0.9 , 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
    r_grid = np.array([0.0326, 0.0270, 0.0234, 0.0204, 0.0182, 0.0167, 0.0154, 0.0143, 0.0134, 0.0126, 0.0119, 0.0112, 0.0105, 0.0099, 0.0093, 0.0087, 0.0082, 0.0076, 0.0070, 0.0064, 0.0057, 0.0050, 0.0042])
    logg_grid = np.log10( (c.G*m_grid*u.Msun/(r_grid*u.Rsun)**2).cgs.value ) # g = GM / r^2
    this_m_interp = np.interp(logg, logg_grid, m_grid) # the mass assumed by the model (corresponding to the log g)
    
    wl_wd, flux_wd = dat.T # surface flux, assuming M = m_interp
    M_wd = this_m_interp*u.Msun
    R_wd_model =  np.sqrt(c.G*M_wd/(10**logg*u.cm/u.s**2)).to(u.Rsun) # calculates correponding radius to the mass assumed by model, r = sqrt(GM/g)
    flux_wd_obs = (R_wd_model/(d_pc*u.pc)).cgs.value**2*flux_wd # (r/d)^2 * flux 
    
    # now scale to the radius it should have, given the mass we want
    r_assume = np.interp(mass_assume, m_grid, r_grid) # mass-radius relation
    flux_wd_obs_this_mass = flux_wd_obs*(r_assume/R_wd_model.value)**2 
    return wl_wd,  flux_wd_obs_this_mass

def kurucz_spectrum_at_a_given_distance(path = 'data/kurucz_bstars/kurucz_Teff_6250_logg_1.5.fits', 
    radius_Rsun = 0.4, parallax = 1.5, vsini = 5):
    '''
    # erg/s/cm^2/AA.  
    '''
    from astropy.io import fits
    # from lb1_tlusty import rotBroad, vmacBroad
    
    c = 2.99792458*1e18 # in AA/s
    pc_cm = 3.086e18 
    G_cgs =  6.67408e-08 #cm3 g-1 s-2
    dist_pc = 1000/parallax
    Msun_g = 1.98848e33
    Rsun_cm = 6.957e+10

    if '.fits' in path: # the BOSZ files 
        hdu = fits.open(path)
        model = hdu[1].data
        
        # Eddington Flux Hl at the stellar surface [erg s-1 cm-2 A-1]. Note that different sources have different units. e.g. https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas is the same but without the factor of pi.
        # wl, Flam_surface, cont_surface = hdu[1].data['Wavelength'], hdu[1].data['SpecificIntensity'], hdu[1].data['Continuum']
        wl, Flam_surface, cont_surface = model.Wavelength.byteswap().newbyteorder(), model.SpecificIntensity.byteswap().newbyteorder(), model.Continuum.byteswap().newbyteorder()
        hdu.close()
        # Factor of pi explained in https://archive.stsci.edu/prepds/bosz/
        mult_factor = np.pi*(radius_Rsun*Rsun_cm)**2/(dist_pc*pc_cm)**2  # pi * (R/d)^2
    elif '.dat.txt' in path: # BT Settle files from svo2
        dat = np.genfromtxt(path)
        wl, Flam_surface, cont_surface = dat.T[0], dat.T[1], dat.T[1] # no cont provided
        mult_factor = (radius_Rsun*Rsun_cm)**2/(dist_pc*pc_cm)**2  
        
    else: # the self-calculated Kurucz models
        dat = np.genfromtxt(path)
        wl, Flam_surface, cont_surface = 10*dat.T[0], dat.T[1], dat.T[2]
        Flam_surface, cont_surface = Flam_surface/10, cont_surface/10 # per AA instead of per 
        mult_factor = 4*np.pi*(radius_Rsun*Rsun_cm)**2/(dist_pc*pc_cm)**2  
        
        
    f_lambda, cont = Flam_surface*mult_factor, cont_surface*mult_factor

    return wl, f_lambda, cont