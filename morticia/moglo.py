__author__ = 'DGriffith'

import numpy as np
# import pandas as pd
# import xray
# import scipy.interpolate
from morticia import Q_  # Import pint quantity class

import warnings
#  MORTICIA globals and utility functions


# The following are the long names that will be provided for plotting purposes and netCDF conventions
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
    'pixapx': 'Pixel Aperture in x',
    'pixapy': 'Pixel Aperture in y',
    'pitchy': 'Pixel Pitch in y',
    'asr': 'Absolute Spectral Response',
    'sqe': 'Spectral Quantum Efficiency',
    'chn': 'Spectral Channel Number',
    'rsr': 'Relative Spectral Response',
    'srf': 'Spectral Response Function',
    'phe': 'Photoelectrons',
    'dn': 'Digital Numbers',
    'wellcap': 'Pixel Well Capacity',
    'readnoise': 'RMS Readout Noise',
    'darkcurr': 'Pixel Dark Current',
    'psnu': 'Photo-Response Non-Uniformity',
    'dsnu': 'Dark Signal Non-Uniformity',
    'tref': 'Reference Temperature',
    'temp': 'Temperature',
    'deltat': 'Temperature Delta',
    'nyqx': 'Nyquist Frequenct in x',
    'nyqy': 'Nyquist Frequency in y',
    'bitdepth': 'Bit Depth',
    'dgain': 'Digital Gain',
    'doffset': 'Digital Offset',
    'dnoise': 'RMS Digital-Equivalent Noise'
}

# The following are intended to be the standard names as per the Climate and Forecast (CF) convention
standard_names = {

}

# Default units are the units that are preferred and to which MORTICIA will convert if other units are given
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
    'fldo': 'deg',
    'phe': 'e',
    'dn': 'count',
    'pitchx': 'mm',
    'pitchy': 'mm',
    'pixapx': 'mm',
    'pixapy': 'mm',
    'wellcap': 'e',
    'readnoise': 'e',
    'darkcurr': 'e/s',  # per pixel
    'prnu': '%',  # Must actually be in percentage  This is a special case and needs special handling
    'dsnu': '%',  # Must actually be in percentage  This is a special case and needs special handling
    'tref': 'degC',
    'temp': 'degC',
    'deltat': 'delta_degC',
    'nyqx': '1/mm',
    'nyqy': '1/mm',
    'bitdepth': 'bit',
    'dgain': 'e/count',
    'doffset': 'count',
    'dnoise': 'count'
}

class Scalar(object):
    """ The Scalar class is for representation of scalar numeric values, together with units of measure, a
    mnemonic, being one of those listed in the MORTICIA long_name vocabulary.

    A feature of the Scalar class is that the data attribute can only be set using a.data = [value, units]

    """

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data_and_units):
        value = Q_(np.asarray(data_and_units[0], dtype=np.float64), data_and_units[1]) # Will blow up if units unknown to pint
        self.__units = data_and_units[1]
        if self.mnemonic in default_units:
            # Convert units
            value = value.to(default_units[self.mnemonic])  # Will blow up if units cannot be converted to default
            self.__units = default_units[self.mnemonic]
        self.__data = value.magnitude

    @property
    def units(self):
        return self.__units

    @units.setter
    def units(self, value):
        warnings.warn('Attempt to alter scalar units ignored.')

    @property
    def mnemonic(self):
        return self.__mnemonic

    @mnemonic.setter
    def mnemonic(self, value):
        warnings.warn('Attempt to alter scalar mnemonic ignored.')

    def __init__(self, mnemonic, data, units, attrs=None):
        """ Construct a scalar value using a MORTICIA mnemonic (from the long_name vocabulary), the scalar value
        and the units of measure.

        :param mnemonic: A string such as 'efl' for Effective Focal Length, which is defined in the MORTICIA
            long_name vocabulary
        :param data: The scalar numeric value to be stored
        :param units: The units string of the data. This string must be known to the Python Pint package.
        :param attrs: A user-defined dictionary of any other attributes for the scalar value. The mnemonic will
            be looked up in moglo.long_name and the long_name will be added as an attribute.
        :return:
        """
        self.__mnemonic = mnemonic
        self.__units = units
        if attrs is None:
            self.attrs = {}
        else:
            self.attrs = attrs
        try:
            self.attrs['long_name'] = long_name[mnemonic]  # Will blow up if the mnemonic is not known
        except KeyError:
            raise KeyError('Unknown scalar mnemonic ' + mnemonic + ' provided in Scalar instantiation.')
        value = Q_(np.asarray(data, dtype=np.float64), units)  # Will blow up if units unknown to pint
        if mnemonic in default_units:
            # Convert units
            value = value.to(default_units[mnemonic])  # Will blow up if units cannot be converted to default
            self.__units = default_units[mnemonic]
        self.__data = value.magnitude

    def __repr__(self):
        return self.attrs['long_name'] + ' : ' + str(self.__data) + ' ' + self.__units


# Define a global exception for unit mismatch
class UnitMismatch(Exception):
    pass

# Define a global exception for missing units
class MissingUnits(Exception):
    pass

class MissingDataArrayAxis(Exception):
    pass