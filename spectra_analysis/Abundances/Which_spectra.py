#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 13:38:12 2022

@author: natsukoyamaguchi
"""
import matplotlib.pyplot as plt
from os import listdir
from scipy import optimize
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian1DKernel, convolve
from PyAstronomy import pyasl

norm = 'rolling'
window = 300 # in AA
col = 'r' # just determines which slice

