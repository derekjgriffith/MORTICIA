# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Perform the standard numpy and units imports
import numpy as np
import matplotlib.pyplot as plt
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity
import matplotlib.pyplot as plt
%matplotlib inline

# <codecell>

import sys
sys.path.append('..')  # Remove this call if the optics module is already installed elsewhere 
#import optics
# This notebook is used for development/testing of the Optics module, so auto reload the Optics module if it changes
%load_ext autoreload
%aimport rad
%autoreload 1

# <codecell>

# Illustrate filter generation using rad.srfgen
w,y,wn,wu = rad.srfgen(center=500, fwhm=10, shape='gauss') # Center wavelength and full width at half maximum default to nm
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Gaussian filter on log scale with some out-of-band leackage between 400 and 600 nm
w,y,wn,wu = rad.srfgen(center=500, fwhm=10, shape='gauss', yedge=1e-7, oob=1e-7, wvmin=400, wvmax=600)
plt.semilogy(wu, y)  # plot against wavelength in microns
plt.xlabel('Wavelength [micron]')

# <codecell>

w,y,wn,wu = rad.srfgen(500, 10, shape='bartlett', wvmin=480, wvmax=520)
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Bartlett (triangular) again with insertion of center flat region of width 10 nm
w,y,wn,wu = rad.srfgen(500, 10, shape='bartlett', centerflat=10, wvmin=480, wvmax=520)
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Generate in wavelength space (nm by default) and plot in wavenumber space
w,y,wn,wu = rad.srfgen(500, 100, shape='welch') # center and full width at half max default to nm
plt.plot(wn, y)  #  plot against wavenumber per cm, note how the filter is obviously skewed in wavenumber space if 
# sufficienty wide
plt.xlabel('Wavenumber [cm^-1]')

# <codecell>

# Similar, but this time, generate the filter in wavenumber scale and also plot in wavenumber scale
w,y,wn,wu = rad.srfgen(20000.0, 2000.0, shape='welch', units='cm^-1')
plt.plot(wn, y)  #  plot against wavenumber per cm
plt.xlabel('Wavenumber [per cm]')

# <codecell>

w,y,wn,wu = rad.srfgen(500, 10, shape='cosine')
plt.plot(w, y)
plt.xlabel('Wavenumber [nm]')

# <codecell>

w,y,wn,wu = rad.srfgen(500, 10, shape='cos^2', wvmin=480, wvmax=520, oob=0.001)
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Same, plotted on a log scale
plt.semilogy(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Same as previous, but open a central flat region of 10 nm width
# Note that the full width at half max is now 20 nm (fwhm + centerflat)
w,y,wn,wu = rad.srfgen(500, 10, shape='cos^2', centerflat=10, wvmin=480, wvmax=520, oob=0.001)
plt.semilogy(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

w,y,wn,wu = rad.srfgen(605, 10, shape='tophat', wvmin=580, wvmax=640)
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Opening a centre flat region also works with tophats (box) although you could just increase fwhm
w,y,wn,wu = rad.srfgen(605, 10, shape='tophat', centerflat=20, wvmin=580, wvmax=640)
plt.plot(w, y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Show the tophat function, which defines a tophat using only a few points
w,y,wn,wu = rad.tophat(605, 10, delta=0.001, wvmin=350, wvmax=750)
plt.plot(w,y)
plt.xlabel('Wavelength [nm]')

# <codecell>

# Create 3 MODTRAN-style flt filters/SRFs with different postions, widths and shapes
filt = rad.Flt('My Special Filters', filterheaders = ['a', 'b', 'c'], centers = [500, 600, 700], fwhms = [10, 20, 30],
               shapes=['gauss', 'cos^2', 'welch'])

# <codecell>

# Write the filters/SRFs in MODTRAN .flt format to a text file
filt.write('SpecialFilters')

# <codecell>

# Plot the filters
filt.plot()

# <codecell>

# Create some filters with just a few points
filt2 = rad.Flt('Hand Filters', filterheaders=['a', 'b'], filters=[np.array([[300, 0.5],[400, 1.0]]),
                                                                   np.array([[300, 0.3],[400, 0.2]])])

# <codecell>

filt2.plot()

# <codecell>

filt2

# <codecell>

CIE = rad.Flt('CIEXYZ2.flt')

# <codecell>

CIE.plot()

# <codecell>

CIE

# <codecell>

Aviris = rad.Flt('')

# <codecell>


# <codecell>

Aviris.write('NewAvirisData.flt')

# <codecell>

import numpy as np
from astropy.nddata import NDData

# <codecell>

a = np.array([1.1,2.2,3.3])

# <codecell>

b = NDData(a)

# <codecell>

c = NDData([5.5,6.6,7.7])

# <codecell>

c

# <codecell>

b

# <codecell>

b+c

# <codecell>

import numpy as np
import pandas as pd
import xray

# <codecell>


