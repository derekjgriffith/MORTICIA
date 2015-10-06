# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Optics Tutorials

# <markdowncell>

# The following tutorials show the usage of various lens-related functions in the optics.py module, which implements a range of functions related to optical systems and the human eye. Many of these functions are used in MORTICA.

# <codecell>

# Perform the standard numpy and units imports
import numpy as np
import matplotlib.pyplot as plt
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity
import matplotlib.pyplot as plt
%matplotlib inline
# This notebook is used for development/testing of the Optics module, so auto reload the Optics module if it changes
#%load_ext autoreload
#%aimport optics
#%autoreload 1
import sys
sys.path.append('../sensor')  # Remove this call if the optics module is already installed elsewhere 
import optics

# If you are pacticipating in development of this notebook, please start the notebook with the --script option

# <codecell>

SpatialFrequencies = Q_(np.arange(0.0, 2200.0, 50.0), '1/millimeter') # cycles per millimeter
Wavelengths = Q_([0.001, 0.0005], 'millimeter') # Wavelengths expressed in mm as well to get valid MTF calculation
FocalRatios = Q_(1.0, '')
MTF = optics.mtf(spf=SpatialFrequencies, wvl=Wavelengths, fno=FocalRatios)

# <codecell>

plt.figure(figsize=(10,8))
plt.rc('text', usetex=True)  # Use TeX to render the labels in the plot
plt.rc('font', family='serif', size=15)  # Default to a serif font 
plt.plot(SpatialFrequencies, MTF)
plt.xlabel('Spatial Frequency $[{:~L}]$'.format(SpatialFrequencies.units))
plt.ylabel('MTF')
plt.title('Modulation Transfer Function')
plt.grid()

# <codecell>

SpatialFrequencies = Q_(np.arange(0.0, 1010.0, 10.0), '1/millimeter') # cycles per millimeter
Wavelengths = Q_(0.001, 'millimeter') # Wavelengths expressed in mm as well to get valid MTF calculation
FocalRatios = [1.0, 2.0]
Obscuration = 0.4  # Aperture is circular with a centred circular obscuration of 40% of the diameter
MTF = optics.mtf(spf=SpatialFrequencies, wvl=Wavelengths, fno=FocalRatios)
MTF = optics.mtf_obs(spf=SpatialFrequencies, wvl=Wavelengths, fno=FocalRatios, obs=Obscuration)

# <codecell>

plt.figure(figsize=(10,8))
plt.rc('text', usetex=True)  # Use TeX to render the labels in the plot
plt.rc('font', family='serif', size=15)  # Default to a serif font
plt.plot(SpatialFrequencies, MTF)
plt.xlabel('Spatial Frequency $[{:~L}]$'.format(SpatialFrequencies.units))
plt.ylabel('MTF')
plt.title('Modulation Transfer Function, Obscuration ' + str(Obscuration*100) + ' \%')
plt.legend(['F/' + str(fno) for fno in FocalRatios])
plt.grid()

# <codecell>


