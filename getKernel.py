# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 19:43:01 2016

@author: luofeng
"""
from __future__ import division
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal

mu = 0
oldR = 30000
newR = 1800
sigma = ((5000./newR)**2-(5000./oldR)**2)**0.5
print 'sigma =', sigma
x = np.linspace(-20, 20, 41)
print x
pi = math.pi
y = []


def GetGaussianKernel(x, mu, sigma):
#    f = math.exp(-(x-mu)**2./(2.0*(sigma**2)))
    f = (1./(sigma*((2.*pi)**0.5)))*math.exp(-(x-mu)**2./(2.0*(sigma**2)))
    return f

for a in x:
    y.append(GetGaussianKernel(a, mu, sigma))

# y_uni=[l/sum(y) for l in y]
# print y
plt.figure(1)
plt.plot(x, y, '-o')
plt.show()

dir = './spec/PHOENIX-ACES-AGSS-COND-2011_Z-0.0/'
filename = 'lte06000-1.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
hl = fits.open(dir+filename)

spec = Table(hl[1].data)
wave = spec['wave']
flux = spec['flux']
kernel = y
flux1 = signal.convolve(flux, kernel, 'same')
plt.figure(2)
# plt.subplot(211)
plt.plot(wave, flux, 'b-')
# plt.subplot(212)
plt.plot(wave, flux1, 'r-')
plt.show()
