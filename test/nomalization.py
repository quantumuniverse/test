# -*- coding: utf-8 -*-

from __future__ import division
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os


def myConv(p, windowsize):
    oldR = 30000
    newR = 1000
    sigma = np.sqrt(pow((5000/newR), 2)-pow((5000/oldR), 2))
    gaussian = np.linspace(0, windowsize, windowsize+1)  # 2.718281828459045
    gaussianF = np.power(np.e, -np.power(gaussian-windowsize/2, 2)/(2*sigma*sigma))
    gaussianF = gaussianF/np.sum(gaussianF)
    postFlux = signal.convolve(p, gaussianF, 'same')
    return postFlux


def mask0(wave, fluxa, fluxfit):
    deli = []
    for i in xrange(len(fluxfit)):
        if np.round(fluxfit[i], 9) == 0.000000000:
            deli.append(i)
    return np.delete(wave, deli), np.delete(fluxa, deli), np.delete(fluxfit, deli)


def normalization(wave, flux, size=3):
#    fluxfilter = signal.medfilt(flux, kernel_size = size)
    fluxfilter = myConv(flux, size)
    wave, fluxa, fluxfit = mask0(wave, flux, fluxfilter)
    divarr = fluxa / fluxfit
    return wave, divarr


dir = '/media/Data/Fits/PHOENIX_A1FITS/PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0'
filename = '/lte06000-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
hl = fits.open(dir+filename)
flux_ph = hl[0].data
naxis1 = hl[0].header['NAXIS1']
crval1 = hl[0].header['CRVAL1']
cdelt1 = hl[0].header['CDELT1']
wave_ph = np.arange(naxis1)*cdelt1 + crval1
indices_ph = np.where((wave_ph >= 3900) & (wave_ph <= 8800))
wave_ph1, flux_ph1 = wave_ph[indices_ph], flux_ph[indices_ph]
wave_ph2, flux_ph2 = normalization(wave_ph1, flux_ph1, 10)


rootdir = "/media/Data/Fits/LAMOST_DR4_snrg_gt_20_lt_100"
for parent, dirnames, filenames in os.walk(rootdir):
    for filename in filenames:
        filename_full = os.path.join(parent, filename)
        hl1 = fits.open(filename_full)
        spec_la = hl1[0].data
        wave_la = spec_la[2]
        flux_la = spec_la[0]
        indices_la = np.where((wave_la >= 3900) & (wave_la <= 8800))
        wave_la1, flux_la1 = wave_la[indices_la], flux_la[indices_la]
        wave_la2, flux_la2 = normalization(wave_la1, flux_la1, 10)

        plt.figure(figsize=(15, 20), dpi=120)
        p1 = plt.subplot(511)
        #p1.axis([3600, 9000, 0, 1.4])
        p2 = plt.subplot(512)
        #p1.axis([3600, 9000, 0, 2])
        p1.plot(wave_ph1, flux_ph1, 'r-')
        p2.plot(wave_ph2, flux_ph2, 'y-')

        p3 = plt.subplot(513)
        p4 = plt.subplot(514)
        p3.plot(wave_la1, flux_la1, 'b-')
        p4.plot(wave_la2, flux_la2, 'g-')
        corr = signal.correlate(flux_la1, flux_ph1, mode='same')
        p5 = plt.subplot(515)
        p5.plot(corr, 'g-')
        plt.savefig('./pics/'+os.path.splitext(filename)[0]+'.png')
        plt.close()


'''
def mask0(wave, fluxa, fluxfit):
    deli = []
    for i in xrange(len(fluxfit)):
        if np.round(fluxfit[i], 9) == 0.000000000:
            deli.append(i)
    return np.delete(wave, deli), np.delete(fluxa, deli), np.delete(fluxfit, deli)


def medfilt_contiu(wave, flux, size=3):
    fluxfilter = signal.medfilt(flux, kernel_size = size)
    wave, fluxa, fluxfit = mask0(wave, flux, fluxfilter)
    divarr = fluxa / fluxfit
    return wave, divarr
'''
