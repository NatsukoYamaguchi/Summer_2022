#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 13:33:10 2022

@author: natsukoyamaguchi
"""
from Which_spectra import *

# folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances'
save_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances' + '/plots' + '/' + norm + '/' + str(window) +  '/' + col
Path(save_dir).mkdir(parents=True, exist_ok=True)

zoom = 'y'
element = 'O'

if element == 'O':
    abundances = [-1.5, -0.5, 0.0, 0.1]
elif element == 'Na' or element == 'N':
    abundances = [-0.2, 0.0, 0.5, 1.0]
elif element == 'C':
    abundances = [-1.5, -0.5, -0.1, 0.2]

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances/spectra_norm' 
if norm == 'poly':
    files = sorted(listdir(folder_dir + '/' + norm + '/' + col)) # alphabetical order
elif norm == 'rolling':
    files = sorted(listdir(folder_dir + '/' + norm + '/' + str(window) + '/' + col)) # alphabetical order

files_plot = []
for i in range(len(abundances)):
    for j in range(len(files)):
        if files[j] == 'at12_t04700g4.80_solar_' + element + str(abundances[i]) + '.spec.gz.csv.csv':
            files_plot.append(files[j])
            
color = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
for i in range(len(files_plot)):
    label = '[' + element + '/H] = ' + str(abundances[i])
    dat = pd.read_csv(folder_dir + '/' + norm + '/' + str(window) + '/' + col + '/' + files_plot[i])
    plt.plot(dat['wavelength'], dat['flux'], label = label, lw = 0.3, c = color[i])

plt.ylabel('normalized flux')
plt.xlabel('wavelength (AA)')
plt.legend()
# plt.ylim(0.2, 1.1)

if zoom == 'y':
    if col == 'b':
        plt.xlim(4000,4500)
    elif col == 'r':
        plt.xlim(5600, 6100)
    plt.savefig(save_dir + '/' + element + '_zoom.png', dpi=300)
else:
    plt.savefig(save_dir + '/' + element + '.png', dpi=300)


#%% Comparison to data 
   
from Which_spectra import *
    
save_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances' + '/plots' + '/' + norm + '/' + str(window) +  '/' + col
Path(save_dir).mkdir(parents=True, exist_ok=True)

element = 'Na'
abundances = [-0.2, 0.0, 1.0]
    
times = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/' + col + '_times.csv')
diff_eclipse = abs(0.5 - times['phase']).to_numpy()
idx_eclipse = np.argwhere(diff_eclipse == np.min(diff_eclipse))[0][0]
jdmid_eclipse = times['jdmid'][idx_eclipse]
phase_eclipse = times['phase'][idx_eclipse]

diff = abs(0.25 - times['phase']).to_numpy()
idx = np.argwhere(diff == np.min(diff))[0][0]
jdmid = times['jdmid'][idx]
phase = times['phase'][idx]

shift_list = pd.read_csv('/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/outputs/shift_' + col + '.csv')
shift = np.exp(shift_list['log_shift'][idx]) - 1
shift_eclipse = np.exp(shift_list['log_shift'][idx_eclipse]) - 1

fn = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/rolling/300/' + col + '_' + str(jdmid) + '.csv'
dat = pd.read_csv(fn)

fn_2 = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/spectra_norm/rolling/300/' + col + '_' + str(jdmid_eclipse) + '.csv'
dat_2 = pd.read_csv(fn_2)

lbl = 'dat, phase = {0:.3f}'.format(phase)
lbl_2 = 'dat, phase = {0:.3f}'.format(phase_eclipse) + ' (eclipse)'

plt.plot(dat['wavelength'] - dat['wavelength']*shift, dat['flux_norm'], label = lbl, lw = 0.3, color = '#5E5E5E', linestyle = '-')
plt.plot(dat_2['wavelength'] - dat_2['wavelength']*shift_eclipse, dat_2['flux_norm'], label = lbl_2, lw = 0.3, color = '#5E5E5E', linestyle = '--')

folder_dir = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis/Abundances/spectra_norm' 
files = sorted(listdir(folder_dir + '/' + norm + '/' + str(window) + '/' + col)) # alphabetical order
files_plot = []

for i in range(len(abundances)):
    for j in range(len(files)):
        if files[j] == 'at12_t04700g4.80_solar_' + element + str(abundances[i]) + '.spec.gz.csv.csv':
            files_plot.append(files[j])

color = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
for i in range(len(files_plot)):
    label = '[' + element + '/H] = ' + str(abundances[i])
    dat = pd.read_csv(folder_dir + '/' + norm + '/' + str(window) + '/' + col + '/' + files_plot[i])
    plt.plot(dat['wavelength'], dat['flux'], label = label, lw = 0.3, c = color[i])

plt.ylabel('normalized flux')
plt.xlabel('wavelength (AA)')
plt.legend()

plt.xlim(5825, 5950)
# plt.xlim(8175, 8210)
# plt.ylim(0.6, 1.2)
plt.savefig(save_dir + '/Na_5900.png', dpi = 300, bbox_inches = "tight")



