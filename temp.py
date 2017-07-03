# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
'''
#%%
import os
os.chdir('/media/luofeng/00060DB10005A1A0/Data/PHOENIX-ACES-AGSS-COND-2011_Z-0.0')


from astropy.io import fits
hl = fits.open('./lte02500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

from astropy.table import Table
spec = Table(hl[1].data)

import matplotlib.pyplot as plt
plt.plot(spec['wave'],spec['flux'])

os.chdir('/home/luofeng/spyderworkspace/test')
plt.savefig("spectrum.pdf");
'''

#%%
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

'''
os.chdir('/media/luofeng/00060DB10005A1A0/Data/PHOENIX-ACES-AGSS-COND-2011_Z-0.0')
hl = fits.open('./lte02500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')
spec = Table(hl[1].data)

x = spec['wave']
y = spec['flux']

fig = plt.figure();
ax = fig.add_subplot(1,1,1)
ax.plot(x,y)

fig.tight_layout()
fig.savefig("spectrum.png");
'''

def GetFileList(dir, fileList):
    newDir = dir
    if os.path.isfile(dir):
        fileList.append(dir.decode('gbk'))
    elif os.path.isdir(dir):  
        for s in os.listdir(dir):
            #如果需要忽略某些文件夹，使用以下代码
            #if s == "xxx":
                #continue
            newDir=os.path.join(dir,s)
            GetFileList(newDir, fileList)  
    return fileList
# /media/luofeng/00060DB10005A1A0/Data/LH_ELODIE_31_23/H				
fileList = GetFileList('/media/D/Data/LL_ELODIE_31_23/L', [])

matplotlib.use('Agg')
os.chdir('/media/D/Data/LL_ELODIE_31_23_spec_pics/L')
for f in fileList:
	fileName = os.path.split(f)[1];
	print os.path.splitext(fileName)[0];
	hl = fits.open(f)
	naxis1 = hl[0].header['NAXIS1']
	crval1 = hl[0].header['CRVAL1']
	crpix1 = hl[0].header['CRPIX1']
	cdelt1 = hl[0].header['CDELT1']
	wave = crval1+(np.arange(naxis1)+1-crpix1)*cdelt1
	flux = hl[0].data
	plt.plot(wave,flux)
	plt.savefig("spectrum_"+ os.path.splitext(fileName)[0] +".png")
	plt.close()
#fig.add_subplot(10,1,1)

#fig.show()
