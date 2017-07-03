# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path


def spec_superposition(rootdir, childdir):
    flux = np.zeros(45001)
    wave_ph = np.zeros(45001)
    for parent, dirnames, filenames in os.walk(rootdir + childdir):
        for filename in filenames:
            filename_full = os.path.join(parent, filename)
            hl = fits.open(filename_full)
            flux_ph = hl[0].data
            naxis1 = hl[0].header['NAXIS1']
            crval1 = hl[0].header['CRVAL1']
            cdelt1 = hl[0].header['CDELT1']
            wave_ph_ori = np.arange(naxis1)*cdelt1 + crval1
            wave_ph = wave_ph_ori
            indices_ph = np.where((wave_ph_ori >= 3900) & (wave_ph_ori <= 8400))
            wave_ph, flux_ph = wave_ph_ori[indices_ph], flux_ph[indices_ph]
            flux = flux + flux_ph
    return wave_ph, flux


def read_spec_ph(filename_full):
    hl = fits.open(filename_full)
    flux_ph = hl[0].data
    naxis1 = hl[0].header['NAXIS1']
    crval1 = hl[0].header['CRVAL1']
    cdelt1 = hl[0].header['CDELT1']
    wave_ph_ori = np.arange(naxis1)*cdelt1 + crval1
    wave_ph = wave_ph_ori
    indices_ph = np.where((wave_ph_ori >= 3900) & (wave_ph_ori <= 8400))
    wave_ph, flux_ph = wave_ph_ori[indices_ph], flux_ph[indices_ph]
    return wave_ph, flux_ph


def read_spec_la(filename_full):
    hl = fits.open(filename_full)
    spec_la = hl[0].data
    wave_la = spec_la[2]
    flux_la = spec_la[0]
    indices_ph = np.where((wave_la >= 3900) & (wave_la <= 8400))
    wave_la, flux_la = wave_la[indices_ph], flux_la[indices_ph]
    return wave_la, flux_la


def traverse_spec(filename, filename_full):
    wave, flux = read_spec_la(filename_full)
    indices_blue = np.where((wave >= 3900) & (wave <= 5800))
    indices_red = np.where((wave > 5800) & (wave <= 8400))
    wave_blue, wave_red = wave[indices_blue], wave[indices_red]
    flux_blue, flux_red = flux[indices_blue], flux[indices_red]
    
    plt.figure(figsize=(24, 12), dpi=70)
    p1 = plt.subplot(311)
    plt.plot(wave, flux, 'b-')
    
    c_blue = np.polyfit(wave_blue, flux_blue, 6)
    c_red = np.polyfit(wave_red, flux_red, 4)
    flux_poly_blue = np.polyval(c_blue, wave_blue)
    flux_poly_red = np.polyval(c_red, wave_red)
    
    plt.plot(wave_blue, flux_poly_blue, 'r-')
    plt.plot(wave_red, flux_poly_red, 'y-')
    plt.xlim(3700, 8600)
    plt.grid(True)
    
    p2 = plt.subplot(312, sharex=p1)
    flux_poly = np.append(flux_poly_blue, flux_poly_red)
    flux_trav = flux / flux_poly
    plt.xlim(3700, 8600)
    plt.grid(True)
    plt.plot(wave, flux_trav, 'g')
    #plt.plot(wave[points], flux_trav[points], 'c.')
    
    p3 = plt.subplot(313)
    nums_bin = np.arange(0.0, 2.0, 0.01)
    n, bins, patches = p3.hist(flux_trav, bins=nums_bin, color='y')
    plt.grid(True)
    # plt.tight_layout()
    m = 0
    for i in np.arange(95, 106, 1):
        m = m + n[i]
    x_label_p3 = "m / n = " + str(round(float(m)/len(flux_trav), 2))
    p3.set_xlabel(x_label_p3)  
    
    #plt.savefig('./spline_pics_lamost/'+os.path.splitext(filename)[0]+'.png')
    print (filename)
    plt.show()
    #plt.close()


def batch_processing(rootdir, childdir):
    for parent, dirnames, filenames in os.walk(rootdir + childdir):
        for filename in filenames:
            filename_full = os.path.join(parent, filename)
            traverse_spec(filename, filename_full)
            break;
            
dirname = '10000-12000/'#'10000-12000'
rootdir = '/media/Data/Fits/'
#childdir = 'PHOENIX_A1FITS/PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0/' + dirname
childdir = 'LAMOST_DR4_snri_snrg_gt_50/'
filename = 'spec-56749-HD151245N190542V01_sp05-106.fits'
#filename = 'lte10000-2.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
# 10000-12000
# 8000-9800
# 6000-7800
# 4000-5900
# 2300-3900
#/media/Data/Fits/LAMOST_DR4_snri_snrg_gt_50
#
#batch_processing(rootdir, childdir)
filename_full = rootdir + childdir + filename
traverse_spec(filename, filename_full)

'''
SDSS young sun lee, sspp, 09 

spl = UnivariateSpline(wave, flux)
#spl.set_smoothing_factor(1)
xs = np.linspace(min(wave), max(wave), 100)
plt.plot(xs, spl(xs), 'r-', lw=1)
'''
