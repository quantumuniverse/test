# -*- coding: utf-8 -*-

from __future__ import division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os


c1 = 37413.82
c2 = 14387.86
count = 0
rootdir = '/media/Data/Fits/PHOENIX_A1FITS/'
childdir = 'PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0/'

for parent, dirnames, filenames in os.walk(rootdir + childdir):
    for filename in filenames:
        filename_full = os.path.join(parent, filename)
        hl = fits.open(filename_full)
        flux_ph = hl[0].data
        naxis1 = hl[0].header['NAXIS1']
        crval1 = hl[0].header['CRVAL1']
        cdelt1 = hl[0].header['CDELT1']
        wave_ph = np.arange(naxis1)*cdelt1 + crval1
        indices_ph = np.where((wave_ph >= 3900) & (wave_ph <= 8800))
        wave_ph1, flux_ph1 = wave_ph[indices_ph], flux_ph[indices_ph]

        plt.figure(figsize=(15, 9), dpi=98)

        p1 = plt.subplot(211)
        flux_ph2 = (flux_ph1-min(flux_ph1))/(max(flux_ph1)-min(flux_ph1))
        p1.plot(wave_ph1, flux_ph2, 'b-')

        T = int(filename[3:8])
        flux_bds = c1/(((wave_ph1/10000.)**5)*(np.exp(c2/((wave_ph1/10000.)*T))-1))
        flux_bds1 = (flux_bds-min(flux_bds))/(max(flux_bds)-min(flux_bds))
        #print "flux_ph1:"+str(min(flux_ph1))+"\t"+str(max(flux_ph1))+"\t"+str(np.mean(flux_ph1))
        #print "flux_ph1:"+str(max(flux_ph1)-min(flux_ph1))
        #print "flux_bds1:"+str(min(flux_bds1))+"\t"+str(max(flux_bds1))+"\n"
        #print "flux_bds1:"+str(max(flux_bds1)-min(flux_bds1))
        p1.plot(wave_ph1, flux_bds1, 'r-')
        p1.set_title(os.path.splitext(filename)[0])
        p1.set_xlabel("wave")
        p1.set_ylabel("flux")

        p2 = plt.subplot(212)
        #T = int(filename[3:8])
        p2.set_title("Black Body " + str(T) +"K")
        p2.set_xlabel("wave")
        p2.set_ylabel("flux")
        #flux_bds = c1/(((wave_ph1/10000.)**5)*(np.exp(c2/((wave_ph1/10000.)*T))-1))
        p2.plot(wave_ph1, flux_ph1/flux_bds, 'g-')
        #plt.savefig('./pics/'+os.path.splitext(filename)[0]+'.png')
        plt.show()
        plt.close()
        count += 1
        if count == 2:
            break
