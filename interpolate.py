# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:44:33 2016

@author: luofeng
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import interpolate

'''
x = np.arange(0,10)
#y = np.exp(-x/3.0)
y = np.sin(x)
f = interpolate.interp1d(x, y, kind='quadratic', axis=-2)

xnew = np.arange(0, 9, 0.1)
ynew = f(xnew)   # use interpolation function returned by `interp1d`
plt.plot(x, y,'o', xnew, ynew, '*')
plt.show()
'''

dir = './spec/LAMOST_DR3_snr_gt_300/'
filename = 'spec-55916-B5591606_sp09-047.fits'
hl = fits.open(dir+filename)

spec = hl[0].data
wave = spec[2]
flux = spec[0]

plt.figure(1)
plt.plot(wave, flux, 'o')
plt.show()

f = interpolate.interp1d(wave, flux, kind='linear')
x = wave
y = flux
xnew = np.arange(min(wave), max(wave), 0.5)
ynew = f(xnew)
plt.figure(2)
plt.plot(x, y, 'o', xnew, ynew, '-*')
plt.show()
