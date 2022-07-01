#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:01:31 2022

@author: natsukoyamaguchi
"""
import matplotlib.pyplot as plt
from os import listdir
from scipy import optimize
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d

base_directory = '/Users/natsukoyamaguchi/Desktop/Summer_2022/lris_p274_analysis/spectra_analysis'
col = 'b' # r or b 
P = 0.11414278 #days 
norm = 'rolling' # poly or rolling
window = 500 # if norm = rolling, choose window size 

Path(base_directory + '/spectra_norm').mkdir(parents=True, exist_ok=True)
Path(base_directory + '/outputs').mkdir(parents=True, exist_ok=True)
Path(base_directory + '/plots').mkdir(parents=True, exist_ok=True)