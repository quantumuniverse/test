# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
from scipy import interpolate
from matplotlib.widgets import MultiCursor
from pylab import figure, show


c1 = 37413.82
c2 = 14387.86
count = 0
rootdir = '/media/Data/Fits/PHOENIX_A1FITS/'
childdir = 'PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0/'
filename = 'lte12000-2.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
filename_full = rootdir + childdir + filename
hl = fits.open(filename_full)
flux_ph_ori = hl[0].data
naxis1 = hl[0].header['NAXIS1']
crval1 = hl[0].header['CRVAL1']
cdelt1 = hl[0].header['CDELT1']
wave_ph_ori = np.arange(naxis1)*cdelt1 + crval1
indices_ph = np.where((wave_ph_ori >= 3900) & (wave_ph_ori <= 8800))
wave_ph, flux_ph = wave_ph_ori[indices_ph], flux_ph_ori[indices_ph]
# flux_ph_01sta = (flux_ph-min(flux_ph))/(max(flux_ph)-min(flux_ph))
# flux_bds_01sta = (flux_bds-min(flux_bds))/(max(flux_bds)-min(flux_bds))
# flux_bds = c1/(((wave_ph/10000.)**5)*(np.exp(c2/((wave_ph/10000.)*T))-1))
# Custom parameters
a = 500  # step = 100
b = 1.0  # step = 0.005
c = 0.80  # step = 0.05
d = 1.0  # step = 0.05
T = int(filename[3:8])
logg = filename.split('-')[1]

def traverse_spectrum():
    plt.figure(figsize=(24, 12), dpi=75)
    p1 = plt.subplot(311)
    flux_bds_modify_t = c1/((((wave_ph - a)/10000.)**5)*(np.exp(c2/(((wave_ph - a)/10000.)*(T*b)))-1))
    flux_bds_trav = (flux_bds_modify_t*10**11 * c) - 10**15*d
    plt.plot(wave_ph, flux_ph, 'b-', label="a500")
    plt.plot(wave_ph, flux_bds_trav, 'r-', label="a500")
    xlabel_p1 = "a="+str(a)+", b="+str(b)+", c="+str(c)+", d="+str(d)
    p1.set_xlabel(xlabel_p1, fontsize=14)
    p1.set_ylabel("flux", fontsize=10)
    p1.set_title(os.path.splitext(filename)[0], fontsize=14)
    plt.grid(True)

    flux_trav = flux_ph/flux_bds_trav
    p2 = plt.subplot(312, sharex=p1)  # sharex
    plt.plot(wave_ph, flux_trav, 'g-')
    p2.set_xlabel("wave", fontsize=10)
    p2.set_ylabel("flux", fontsize=10)
    plt.xlim(3700, 9000)
    plt.grid(True)

    p3 = plt.subplot(313)
    num_bins = np.arange(0.0, 2.0, 0.005)
    n, bins, patches = p3.hist(flux_trav, bins=num_bins, color='y')
    p3.set_xlabel("flux", fontsize=10)
    p3.set_ylabel("frequency", fontsize=10)
    plt.grid(True)

# plt.tight_layout()
# plt.show()
    filename_pic = os.path.splitext(filename)[0] + "  " + xlabel_p1
    
    plt.savefig('./pics/' + str(count) + filename_pic + '.png')
    plt.close()
    print (xlabel_p1)

for a in np.arange(0, 300, 50):
    for b in np.arange(0.95, 1.06, 0.02):
        for c in np.arange(0.95, 1.06, 0.02):
            for d in np.arange(0.0, 1.1, 0.2):
                count += 1
                print (count)
                traverse_spectrum()

'''
待完成工作：
1.图上添加标签注释
2.循环4个参数，探索跟Log g的关系
3.循环9个文件
4.保存图片
'''
