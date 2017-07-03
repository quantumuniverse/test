#!/usr/bin/env python
import os
import os.path
from astropy.io import fits
import numpy as np  
import pylab as plt  
from scipy import stats
from scipy.interpolate import splrep, splev
#import matplotlib.pyplot as plt  


def read_spec_ph(filename_full):
    hl = fits.open(filename_full)
    flux_ph = hl[0].data
    naxis1 = hl[0].header['NAXIS1']
    crval1 = hl[0].header['CRVAL1']
    cdelt1 = hl[0].header['CDELT1']
    wave_ph_ori = np.arange(naxis1)*cdelt1 + crval1
    #wave_ph = wave_ph_ori
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


def choose_point1(points, flux):
    p = points
    f = flux
    i = 1
    result = []
    result.append(p[0])
    while(i < len(p)-1):
        if(f[p[i-1]] < f[p[i]] and f[p[i]] > f[p[i+1]]):
            result.append(p[i])
        i = i + 1
    result.append(p[-1]) 
    return np.array(result)


def choose_point2(points, flux):
    p = points
    f = flux
    i = 1
    result = []
    result.append(p[0])
    while(i < len(p)-1):
        if(f[p[i-1]] > f[p[i]] and f[p[i]] < f[p[i+1]]):
            i = i + 1
            continue;
        result.append(p[i])
        i = i + 1
    result.append(p[-1]) 
    return np.array(result)


def reject_pixels_3sigma(points, flux_norm):
    print ('+++ in 3sigma +++')
    mu, sigma = stats.norm.fit(flux_norm)
    top = mu + sigma*1.5
    bottom = mu - sigma*3
    print (mu, sigma, bottom, top)
    #indices = np.where((flux_norm[points] < top) & (flux_norm[points] > bottom))
    indices = np.where(flux_norm[points] < top)
    #flux_norm = flux_norm[indices]
    #print (len(w), len(flux_norm))
    print (10, len(indices[0]))
    print ('--- in 3sigma ---')
    print ()
    return points[indices[0]]

    
def choose_point(f):
    print ('+++ in choose_point +++')
    points = np.arange(0, len(f), 1)
    print (1,len(points))
    points = choose_point1(points, f)
    print (2,len(points))
    points = choose_point1(points, f)
    print (3,len(points))
    points = choose_point2(points, f)
    print (4,len(points))
    points = choose_point2(points, f)
    print (5,len(points))
    points = choose_point2(points, f)
    print (6,len(points))
    points = choose_point2(points, f)
    print (7,len(points))
    points = choose_point2(points, f)
    print (8,len(points))
    points = choose_point2(points, f)
    if (f[points[0]] < f[points[1]]):
        points = points[1:]
        '''
    if (f[points[-1]] < f[points[-2]]):
        points = points[:-1]
        '''
    print (9,len(points))
    print ('--- in choose_point ---')
    print ()
    return points

    
def spline_fit(wave, w, f, smoothing_factor):
    print ('+++ in spline_fit +++')
    print (smoothing_factor)
    spl = splrep(w, f, s=smoothing_factor)
    flux_spl = splev(wave, spl)
    print ('--- in spline_fit ---')
    print ()
    return flux_spl


def flux_normalization(flux, flux_fit):
    return flux / flux_fit


def get_image(wave, flux, flux_spl, points, flux_norm, filename):
    print ('+++ in get_image +++')
    plt.figure(figsize=(24, 12), dpi=75)
    p1 = plt.subplot(311)
    plt.plot(wave, flux, 'b-')
    plt.plot(wave, flux_spl, 'r-')
    #plt.plot(wave[points], flux[points], 'y.')
    plt.xlim(3700, 8600)
    plt.grid(True)
    p2 = plt.subplot(312, sharex=p1)
    #print(np.where(f == 0))
    plt.xlim(3700, 8600)
    plt.grid(True)
    plt.plot(wave, flux_norm, 'g-')
    #plt.plot(wave[points], flux_trav[points], 'r.')
    p3 = plt.subplot(313)
    #mu, sigma = stats.norm.fit(flux_trav)
    #plt.plot(wave, stats.norm.pdf(wave, mu, sigma))
    nums_bin = np.arange(0.0, 2.0, 0.01)
    n, bins, patches = p3.hist(flux_norm, bins=nums_bin, color='y')
    plt.grid(True)
    
    m = 0
    for i in np.arange(95, 106, 1):
        m = m + n[i]
    
    x_label_p3 = "m / n = " + str(round(float(m)/len(flux_norm), 2))
    p3.set_xlabel(x_label_p3)  
    
    file = os.path.splitext(filename)[0] + '.png'
    print (filename)
    plt.tight_layout()
    plt.savefig(path + file)
    #plt.show()
    plt.close()
    print ('--- in get_image ---')
    print ()


def traverse_spec(filename, filename_full):
    print ('+++ in traverse +++')
    wave, flux = read_spec_ph(filename_full)
    flux = (flux-min(flux)) / (max(flux)-min(flux)) + 1
    points = choose_point(flux)
    s = 0.2
    flux_spl = spline_fit(wave, wave[points], flux[points], s)
    flux_norm = flux_normalization(flux, flux_spl)
    #get_image(wave, flux, flux_spl, points, flux_norm, filename)
    print ('points '+str(len(points)))
    points = reject_pixels_3sigma(points, flux_norm)
    print ('points '+str(len(points)))
    s = len(points)/2500
    flux_spl = spline_fit(wave, wave[points], flux[points], s)
    flux_norm = flux_normalization(flux, flux_spl)
    get_image(wave, flux, flux_spl, points, flux_norm, filename)
    print ('--- in traverse ---')
    print ()
    return flux_spl


def batch_processing(rootdir, childdir):
    for parent, dirnames, filenames in os.walk(rootdir + childdir):
        for filename in filenames:
            filename_full = os.path.join(parent, filename)
            traverse_spec(filename, filename_full)
            #break;

rootdir = '/media/Data/Fits/'
childdir = 'PHOENIX_A1FITS/PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z-0.0/'#8000-9800/
#childdir = 'LAMOST_DR4_snri_snrg_gt_50/'
'''
#filename = 'lte10000-2.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
filename = 'spec-56811-KP190339N395439V01_sp09-238.fits'

filename_full = rootdir + childdir + filename
traverse_spec(filename, filename_full)
'''
date = '2017-07-01'
path = '/media/Data/Fits/Pictures/' + date + '/spline_pics_phoenix/'
if (not os.path.exists(path)):
    os.makedirs(path)
print (rootdir+childdir)
batch_processing(rootdir, childdir)

'''
tough spectra
spec-56568-EG023131N032619B01_sp10-146
spec-56568-EG023131N032619F07_sp10-146
spec-56570-EG225232N033308V01_sp08-069
spec-56659-HD104509N022916B01_sp16-057
spec-56751-HD174813N381357V01_sp14-208
spec-56811-HD160847N353259B01_sp15-204
spec-56811-KP190339N395439V01_sp01-121
spec-56811-KP190339N395439V01_sp09-238
'''