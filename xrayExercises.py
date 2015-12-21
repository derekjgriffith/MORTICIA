
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
import xray
import matplotlib.pyplot as plt
import scipy.interpolate
get_ipython().magic(u'matplotlib inline')


# In[2]:

# Build a spectral transmission curve
spec_trans = xray.DataArray([ 0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.0, 0.8, 0.0],                      
                   [('wvl', [550., 600, 650, 700, 750, 800, 850, 950, 1000])], 
                   name='trn',
                   attrs={'trn_units': '1', 'wvl_units': 'nm', 'extrap_val': 0.0})
spatial_frequencies = xray.DataArray(np.linspace(0.0, 1000.0, 11), name='spf', attrs={'spf_units': 'cy/mm'})


# In[3]:

# Coordinate axis names are used to plot the data
spec_trans.plot()
plt.grid()


# In[4]:

spec_trans.coords['wvl']


# In[5]:

spec_trans.dims


# In[6]:

spec_trans.attrs['wvl_units']


# In[7]:

# Create a new wavelength grid
wvl = xray.DataArray(np.linspace(500., 1050, 51), name='wvl', attrs={'wvl_units': 'nm'})


# In[ ]:

wvl


# In[ ]:

spc = scipy.interpolate.interp1d(spec_trans.coords['wvl'].data, spec_trans.data, kind='slinear', bounds_error=False, fill_value=0.0)


# In[ ]:

plt.plot(wvl.data, spc(wvl.data))
plt.grid()


# In[ ]:

spec_trans.data


# In[ ]:

wvl


# In[ ]:

x = scipy.interpolate.RegularGridInterpolator((spec_trans.coords['wvl'].data,), spec_trans.data, bounds_error=False, fill_value=0.0)


# In[ ]:

x(wvl.data)


# In[ ]:

spec = xray.DataArray([ 0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.0, 0.8, 0.0],                      
                   [('refl', [550., 600, 650, 700, 750, 800, 850, 950, 1000])], 
                   name='refl',
                   attrs={'refl_units': '1', 'wvl_units': 'nm', 'extrap_val': 0.0})


# In[ ]:

spec*spec_trans


# In[ ]:

w = xray.DataArray([500.,600.,700.], name='wvl', attrs={'units': 'nm'})


# In[ ]:

f = xray.DataArray([0.3, .4, .2], [(w.name, w.data)], name='filt', attrs={'units_filt': '1', 'units_'+w.name: w.attrs['units']})


# In[ ]:

f


# In[14]:

from . import moglo


# In[1]:

import os


# In[4]:

os.path.abspath(os.path.curdir)


# In[ ]:



