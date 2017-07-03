# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
from scipy.interpolate import splrep, splev

dirname = '10000-12000'
rootdir = '/media/Data/Fits/PHOENIX_A1FITS/'
childdir = 'PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0/' + dirname
# filenames[0] = 'lte12000-2.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
# 10000-12000
# 8000-9800
# 6000-7800
# 4000-5900
# 2300-3900

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
            # flux = [flux[i]+flux_ph[i] for i in range(len(flux_ph))]
            flux = flux + flux_ph
            #print (filename)
    return wave_ph, flux


def read_spec(rootdir, childdir, filename):
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
            # flux = [flux[i]+flux_ph[i] for i in range(len(flux_ph))]
            flux = flux + flux_ph
            #print (filename)
    return wave_ph, flux

def choose_point1(flux):
    f = flux
    i = 1
    result = []
    result.append(0)
    while(i < len(f)-1):
        # print 'i=' + str(i)
        if(f[i-1] < f[i] and f[i] > f[i+1]):
            result.append(i)
        i = i + 1
    result.append(len(f)-1) 
    return result


def choose_point2(points, flux):
    p = points
    f = flux
    i = 1
    result = []
    result.append(p[0])
    while(i < len(p)-1):
        # print 'i=' + str(i)
        
        if(f[p[i]] < f[p[i-1]] and f[p[i]] < f[p[i+1]]):
            #print (i)
            i = i + 1
            continue;
        result.append(p[i])
        #print ('i='+str(i)+' p[i]='+str(p[i]))
        i = i + 1
    result.append(p[-1]) 
    return result


def choose_point3(points, flux):
    p = points
    f = flux
    i = 1
    result = []
    result.append(p[0])
    while(i < len(p)-1):
        # print 'i=' + str(i)
        
        if(f[p[i]] < f[p[i-1]] and f[p[i]] < f[p[i+1]]):
            #print (i)
            i = i + 1
            continue;
        result.append(p[i])
        #print ('i='+str(i)+' p[i]='+str(p[i]))
        i = i + 1
    result.append(p[-1]) 
    return result

#def traverse_spec(wave, flux):
wave, flux = spec_superposition(rootdir, childdir)
plt.figure(figsize=(24, 12), dpi=75)
p1 = plt.subplot(311)
plt.plot(wave, flux, 'b-')

#spl = UnivariateSpline(wave, flux)
#spl.set_smoothing_factor(1)

#points = np.linspace(min(wave), max(wave), 10)
points = choose_point1(flux)
print (len(points))
points = choose_point2(points, flux)
print (len(points))
#x = wave[points]
#print len(x)
spl = splrep(wave[points], flux[points], s=0)
flux_spl = splev(wave, spl)
#flux_spl = spl(x)
plt.xlim(3700, 8600)
plt.grid(True)

plt.plot(wave, flux_spl, 'r-')

'''
points = choose_point1(flux_spl)
print (len(points))
points = choose_point2(points, flux_spl)
print (len(points))
#x = wave[points]
#print len(x)
spl = splrep(wave[points], flux_spl[points], s=0)
flux_spl = splev(wave, spl)
#flux_spl = spl(x)
plt.plot(wave, flux_spl, 'c-')
'''
p2 = plt.subplot(312, sharex=p1)
flux_trav = flux / flux_spl
#points = choose_point3(points, flux_trav)
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
    #print (i, m)

'''
#print (round(float(m)/len(flux_trav), 2))
#print (format(float(m)/len(flux_trav), '.2f'))
#print ('%.2f' %(float(m)/len(flux_trav)))
'''

x_label_p3 = "m / n = " + str(round(float(m)/len(flux_trav), 2))
p3.set_xlabel(x_label_p3)  

#for p in patches:
#    print(p)

plt.savefig('./spline_pics/'+dirname+'K01.png')
plt.show()
# plt.close()

'''
spl = UnivariateSpline(wave, flux)
#spl.set_smoothing_factor(1)
xs = np.linspace(min(wave), max(wave), 100)
plt.plot(xs, spl(xs), 'r-', lw=1)
'''