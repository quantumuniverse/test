# -*- coding: utf-8 -*-

from __future__ import division
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

dir = '/media/D/Data/LAMOST_DR4_snrg_gt_20_lt_100'
hl1 = fits.open(dir+'/spec-55862-B6202_sp09-035.fits')
spec_la = hl1[0].data
wave_la = spec_la[2]
flux_la = spec_la[0]
indices_la = np.where((wave_la >= 3900) & (wave_la <= 8800))
wave_la1 = wave_la[indices_la]
flux_la1 = flux_la[indices_la]
plt.figure(figsize=(15,12),dpi=75)
p1=plt.subplot(311)
p1.plot(wave_la1, flux_la1, 'b-')
#plt.savefig('la.png')

dir = '/media/D/Data/PHOENIX-ACES-AGSS-COND-2011_Z-0.0'
filename = '/lte06000-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
hl = fits.open(dir+filename)
spec_ph = Table(hl[1].data)
wave_ph = spec_ph['wave']
flux_ph = spec_ph['flux']
indices_ph = np.where((wave_ph >= 3900) & (wave_ph <= 8800))
wave_ph1 = wave_ph[indices_ph]
flux_ph1 = flux_ph[indices_ph]
p2=plt.subplot(312)
p2.plot(wave_ph1, flux_ph1, 'r-')
#plt.savefig('ph.png')

corr = signal.correlate(flux_la1, flux_ph1, mode='same')
p3 = plt.subplot(313)
p3.plot(corr, 'g-')
plt.savefig('dpi_75_12.png')

'''
print len(wave_la[indices_la]), len(wave_la1)
print len(wave_ph[indices_ph]), len(wave_ph1)
print min(wave_la), max(wave_la)
print min(wave_ph), max(wave_ph)
print 'corr=', len(corr)
print '-------------'
print indices_la[0]
print indices_ph[0]
'''
