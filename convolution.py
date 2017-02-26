# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt


def myConv(p, windowsize):
    oldR = 30000
    newR = 1000
    sigma = np.sqrt(pow((5000/newR), 2)-pow((5000/oldR), 2))
    print sigma
    gaussian = np.linspace(0, windowsize, windowsize+1)  # 2.718281828459045
    gaussianF = np.power(np.e, -np.power(gaussian-windowsize/2, 2)/(2*sigma*sigma))
    gaussianF = gaussianF/np.sum(gaussianF)
    postFlux = scipy.signal.convolve(p, gaussianF, 'same')
#    plt.figure(1)
#    plt.plot(gaussian,gaussianF,'-o')
#    plt.show()
    return postFlux

dir = './spec/PHOENIX-ACES-AGSS-COND-2011_Z-0.0/'
filename = 'lte02300-0.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

hl = fits.open(dir+filename)
spec = Table(hl[1].data)

wave = spec['wave']
flux = spec['flux']

plt.figure(2)
plt1 = plt.plot(wave, flux, '-', label='origin')

flux1 = myConv(flux, 10)
plt2 = plt.plot(wave, flux1, '-', label='w = 10')

flux1 = myConv(flux, 50)
plt3 = plt.plot(wave, flux1, '-', label='w = 50')
plt.xlabel('wavelength')
plt.ylabel('flux')
plt.legend()
plt.show()
