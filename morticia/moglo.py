__author__ = 'DGriffith'

# import numpy as np
# import pandas as pd
# import xray
# import scipy.interpolate

# MORTICIA globals and utility functions


# The following are the long names that will be provided
long_names = {
    'wvl': 'Wavelength',
    'wvn': 'Wavenumber',
    'trn': 'Transmission',
    'rad': 'Radiance',
    'irr': 'Irradiance',
    'fldx': 'Field Position in x',
    'fldy': 'Field Position in y',
    'fldz': 'Defocus',
    'spf': 'Spatial Frequency',
    'mtf': 'Modulation Transfer Function',
    'otf': 'Optical Transfer Function',
    'ptf': 'Phase Transfer Function',
    'fno': 'Focal Ratio',
    'efl': 'Effective Focal Length',
    'lum': 'Luminance',
    'ill': 'Illuminance',
    'flux': 'Optical Flux',
    'flo': 'Field Orientation',  # Typically horizontal/vertical, across/along track, sagittal/tangential
    'obs': 'Obscuration Ratio',
    'pitchx': 'Pixel Pitch in x',
    'pitchy': 'Pixel Pitch in y',
    'asr': 'Absolute Spectral Response',
    'sqe': 'Spectral Quantum Efficiency'
}

# The following are intended to be the standard names as per the Climate and Forecast (CF) convention
standard_names = {

}