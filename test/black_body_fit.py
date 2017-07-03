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
        wave_ph_ori = np.arange(naxis1)*cdelt1 + crval1
        indices_ph = np.where((wave_ph_ori >= 3900) & (wave_ph_ori <= 8800))
        wave_ph, flux_ph = wave_ph_ori[indices_ph], flux_ph[indices_ph]

        plt.figure(figsize=(20, 12), dpi=75)

        p1 = plt.subplot(311)
        flux_ph_01sta = (flux_ph-min(flux_ph))/(max(flux_ph)-min(flux_ph))
        p1.plot(wave_ph, flux_ph_01sta, 'b-')

        T = int(filename[3:8])
        flux_bds = c1/(((wave_ph/10000.)**5)*(np.exp(c2/((wave_ph/10000.)*T))-1))
        flux_bds_01sta = (flux_bds-min(flux_bds))/(max(flux_bds)-min(flux_bds))

        flux_trav = flux_ph/flux_bds
        #p1.set_title(os.path.splitext(filename)[0])
        p1.set_xlabel("wave")
        p1.set_ylabel("flux")
        p1.axis([3700, 9000, 0, 1.05])
        p1.plot(wave_ph, flux_bds_01sta, 'r-')
        p1.grid(True)

        p2 = plt.subplot(312)
        p2.plot(wave_ph, flux_trav, 'g-')
        #p2.set_title("Black Body Fit " + str(T) + "K")
        p2.set_xlabel("wave")
        p2.set_ylabel("flux")
        plt.xlim(3700,9000)
        p2.grid(True)

        p3 = plt.subplot(313)
        num_bins = np.arange(0.3*10**11, 1.4*10**11, 0.005*10**11)
        n, bins, patches = p3.hist(flux_trav, bins=num_bins)
        p3.grid(True)
        #plt.xlim(0.3, 1.4)
        p3.set_xlabel("flux")
        p3.set_ylabel("frequency")
        # plt.savefig('./pics/'+os.path.splitext(filename)[0]+'.png')
        plt.show()
        # print (n), (bins)
        #plt.close()
        count += 1
        if count == 1:
            break
