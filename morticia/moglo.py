__author__ = 'DGriffith'

# import numpy as np
# import pandas as pd
# import xray
# import scipy.interpolate

# MORTICIA globals and utility functions


# The following are the long names that will be provided for plotting purposes and
long_name = {
    'wvl': 'Wavelength',
    'wvn': 'Wavenumber',
    'trn': 'Transmission',
    'rad': 'Radiance',
    'irr': 'Irradiance',
    'specrad': 'Spectral Radiance',
    'specirr': 'Spectral Irradiance',
    'spf': 'Spatial Frequency',
    'mtf': 'Modulation Transfer Function',
    'otf': 'Optical Transfer Function',
    'ptf': 'Phase Transfer Function',
    'fno': 'Focal Ratio',
    'efl': 'Effective Focal Length',
    'lum': 'Luminance',
    'ill': 'Illuminance',
    'flux': 'Optical Flux',
    'fldy': 'Field Position in x',
    'fldx': 'Field Position in y',
    'fldz': 'Defocus',
    'fldo': 'Field Orientation',  # Typically horizontal/vertical, across/along track, sagittal/tangential
    'obs': 'Obscuration Ratio',
    'pitchx': 'Pixel Pitch in x',
    'pitchy': 'Pixel Pitch in y',
    'asr': 'Absolute Spectral Response',
    'sqe': 'Spectral Quantum Efficiency',
    'chn': 'Spectral Channel Number',
    'rsr': 'Relative Spectral Response',
    'srf': 'Spectral Response Function'
}

# The following are intended to be the standard names as per the Climate and Forecast (CF) convention
standard_names = {

}

default_units = {
    'wvl': 'nm',
    'spf': '1/mm',
    'efl': 'mm',
    'rad': 'W/m^2/sr',
    'irr': 'W/m^2',
    'specrad': 'W/m^2/sr/nm',
    'specirr': 'W/m^/nm',
    'fldx': 'mm',
    'fldy': 'mm',
    'fldz': 'mm',
    'flux': 'W',
    'asr': 'A/W',
    'sqe': '',
    'lum': 'cd/m^2',
    'ill': 'lux',
    'flo': 'deg',
}

# Define a global exception for unit mismatch
class UnitMismatch(Exception):
    pass

# Define a global exception for missing units
class MissingUnits(Exception):
    pass

class MissingDataArrayAxis(Exception):
    pass