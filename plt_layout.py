"""
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import datetime

# my fake data
dates = np.array([datetime.datetime(2000,1,1) + datetime.timedelta(days=i) for i in range(365*5)])
data = np.sin(np.arange(365*5)/365.0*2*np.pi - 0.25*np.pi) + np.random.rand(365*5) /3

# creates fig with 2 subplots
fig = plt.figure(figsize=(10.0, 6.0))
ax = plt.subplot2grid((2,1), (0, 0))
ax2 = plt.subplot2grid((2,1), (1, 0))
## plot dates
ax2.plot_date( dates, data )

# rotates labels 
plt.setp( ax2.xaxis.get_majorticklabels(), rotation=-45 ) 

# shift labels to the right
for tick in ax2.xaxis.get_majorticklabels():
    tick.set_horizontalalignment("right")

plt.tight_layout()
plt.show()
"""

'''
import matplotlib.pyplot as plt

X1 = range(0, 50)

Y1 = [num**2 for num in X1] # y = x^2

X2 = [0, 1]

Y2 = [0, 1]  # y = x

Fig = plt.figure(figsize=(8,4)) # Create a `figure' instance

Ax = Fig.add_subplot(10,1,1) # Create a `axes' instance in the figure

Ax.plot(X1, Y1, X2, Y2) # Create a Line2D instance in the axes

Fig.show()

#Fig.savefig("test.pdf")
'''
'''
#%%
import numpy as np
import matplotlib.pyplot as plt
plt.figure(1) 
plt.figure(2) 
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
x = np.linspace(0, 3, 100)
for i in xrange(3):
    plt.figure(1)
    plt.plot(x, np.exp(i*x/3))
    plt.sca(ax1)
    plt.plot(x, np.sin(i*x))
    plt.sca(ax2)
    plt.plot(x, np.cos(i*x))
plt.show()
#plt.figure(1).savefig("tt.pdf")
'''
#%%
import os
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/media/D/Data/PHOENIX-ACES-AGSS-COND-2011_Z-0.0')
#hl = fits.open('./lte02500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')
#spec = Table(hl[1].data)

plt.figure(1) 
x = np.linspace(0, 3, 100)
plt.figure(1)
plt.plot(x, np.exp(5*x/3))
plt.show()
#plt.figure(1).savefig("tt.pdf")