__author__ = 'DGriffith, ARamkilowan'

"""
.. module:: electro
    :synopsis: Realises electronic components in the imaging chain of a remote sensing or surveillance system.
        - Detectors and focal plane arrays
        - Image processing modules
        - Displays

"""

import numpy as np
import pandas as pd
import xray
from morticia import ureg, Q_, U_  # Import the pint units registry from parent
from morticia.tools.xd import *  # Import additional tools for working with xray.DataArray objects
from morticia.moglo import *  # Import global glossary/vocalbulary

def xd_asr2sqe(asr):
    """ Convert absolute spectral response (ASR) to spectral quantum efficiency (SQE).

        The conversion is performed through the Planck relationship

        .. math::
            E=h\\nu

        where :math:`E` is the photon energy, :math:`h` is the planck constant and :math:`\\nu` is the optical
        frequency.
        In terms of optical wavelength in a vacuum, the planck relation is

        .. math::
            E=\\frac{hc}{\\lambda}

        where :math:`c` is the speed of light and :math:`\\lambda` is the wavelength.

    :param asr: An xray.DataArray object providing the absolute spectral response (ASR). The DataArray must have a single
        axis providing the wavelength points, together with the standard attribute 'units' for the wavelength axis.
    :return: Spectral Quantum Efficiency as an xray.DataArray object. Returned wavelengths will be 'nm' in
        the wavelength ('wvl') axis.

    """
    h = 6.626069311e-34  # planck's constant in Joule seconds
    c = 299792458  # speed of light in metres per second
    e = 1.6021765314e-19  # charge on electron in coulombs
    xd_check_convert_units(asr, 'wvl', 'nm')  # in-place conversion using pint
    xd_check_convert_units(asr, 'asr', 'A/W')  # Preferred units for SQE is A/W
    # Get quantum energy
    E = h * c / asr['wvl'] / 1.0e-9  # Wavelength in nm, convert to metres
    photocurrent = asr / e
    sqe = photocurrent * E
    sqe.attrs['units'] = ''
    sqe.attrs['long_name'] = long_name['sqe']
    sqe['wvl'].attrs['units'] = 'nm'
    sqe['wvl'].attrs['long_name'] = 'Wavelength'
    if any(sqe.values > 1.0):
        warnings.warn('electro.xd_asr2sqe computed unphysical quantum efficiencies exceeding 1.')
    return sqe


def xd_sqe2asr(sqe):
    """ Convert spectral quantum efficiency (SQE) to absolute spectral response.
        This is the reverse conversion of that provided by electro.xd_asr2sqe.

    :param sqe: The spectral quantum efficiency (SQE), provided as a xray.DataArray object in which SQE is
        provided as a function of wavelength. The wavelength axis must provide the 'units' attribute.
    :return: Absolute spectral response (ASR) as an xray.DataArray object, havin a single axis of wavelength
        coordinates ('wvl'). Returned units for the wavelength axis will be 'nm'.

    """
    h = 6.626069311e-34  # planck's constant in Joule seconds
    c = 299792458  # speed of light in metres per second
    e = 1.6021765314e-19  # charge on electron in coulombs
    xd_check_convert_units(sqe, 'wvl', 'nm')  # in-place conversion using pint
    # Get quantum energy
    E = h * c / sqe['wvl'] / 1.0e-9  # Wavelength in nm, convert to metres
    photocurrent = sqe / E
    asr = photocurrent * e
    asr.attrs['units'] = 'A/W'
    asr.attrs['long_name'] = long_name['asr']
    asr['wvl'].attrs['units'] = 'nm'
    asr['wvl'].attrs['long_name'] = 'Wavelength'
    return asr


class FocalPlaneArray(object):
    """ Focal plane array detector. This implementation is typically at the chip level. That is, all or most of the
        information for building an FPA object can be found in the chip-level datasheet. The FPAs in question here
        are usually CCD, CMOS, scientific CMOS or electron-multiplying CCD (EMCCD).

        The FPA class can be combined with an image intensifier tube (IIT) to produce an ICCD device.

        This class does not model Time-Delay and Integration (TDI), which is a dynamic imaging process.

        The class does allow for setting of FPA operating temperature and recalculation of dark current
        based on the dark current doubling delta temperature.

    """
    # Define the temperature property and associated setter method
    # This is required because the operating temperature affects the dark current
    @property
    def temperature(self):
        return self.temperature

    @temperature.setter
    def temperature(self, operating_temperature):
        """ Set the current operating temperature of the FocalPlaneArray. This influences the dark current
            according to the dark current doubling temperature. Result will not be reliable if the dark
            current doubling temperature attribute of the FPA is not set correctly.

        :param operating_temperature:
        :return:
        """
        check_convert_units(operating_temperature, 'degC')
        # Set through self dictionary to avoid recursion
        self.__temperature = operating_temperature[0]
        self.set_dark_current()

    @temperature.getter
    def temperature(self):
        return self.__temperature

    def __init__(self, pitch, aperture, pixels, wellcapacity, readnoise, darkcurrent, dsnu, prnu, sqe=None, asr=None,
                 t_ref=(25.0, 'degC'), darkcurrent_delta_t=(7.0, 'delta_degC'), temperature=(25.0, 'degC'),
                 attrs=None):
        """ Constructor for FocalPlaneAArray objects

        :param pitch: The centre-to-centre spacing of the FPA detector elements. This must be a list where the
            first element is the centre-to-centre spacing in the x-direction (or both x and y), the second element is
            the spacing in the y-direction (if different from the x-direction) and the last element is the units
            in which the pitch is provided (string). The units must be pint-recognizable.
        :param aperture: The effective pixel aperture of the FPA detector elements. This must also be a list, providing
            the x-aperture of the pixel, y-aperture (if different from x) and the units (string).
        :param pixels: A list of number of pixels in the x and y directions respectively
        :param wellcapacity: The well capacity of the pixel in number of electrons
        :param readnoise: The RMS read noise in electrons.
        :param darkcurrent: The dark current as a list, giving nagnitude and units. Units can be electrons per pixel
            per second (e/s), or in an amperage per pixel (A/s) or an amperage per unit area of the FPA (e.g. A/cm^2)
        :param dsnu: The dark signal non uniformity provided as a list giving the magnitude in the first element
            and units in the second element. The units can also be e/s, A or A/area. This is the standard deviation.
            The DSNU can also be specified as 0.0 or None.
        :param prnu: The photo-response non-uniformity. This is also a standard deviation and must be provided in units
            of '%'. That is, the prnu input is a list with the magnitude in the first position and the obligatory
            string '%' in the second position.
        :param sqe: This is a xray.DataArray object providing the Spectral Quantum Efficiency of the FPA. SQE must be
            provided with a single axis of wavelength values. The attribute sqe_units must be provided in the
            DataArray attributes and it must be dimensionless (the literal empty string ''). If any of the values
            exceeds unity, the SQE is assumed to be provided in percent. Alternatively, the asr input can be provided,
            but either the sqe or the asr input must be provided and not both, since they can be converted one from
            the other. SQE outside the spectral domain provided is assumed to be zero.
        :param asr: This is a xray.DataArray object providing the Absolute Spectral Response of the FPA. It can be
            provided as an alternative to the SQE. It must have the single axis of wavelength. Units are equivalent
            to A/W (photoelectron current per unit optical flux).
        :param t_ref: The reference temperature at which the dark current and other parameters are specified.
            To be provided as a [value, 'units'] list. Default is [25.0, 'degC'].
        :param darkcurrent_delta_t: The increase in temperature that causes doubling of the dark current.
            Default is [7.0, 'delta_degC'].
        :param temperature: The operating temperature of the FocalPlaneArray. Default [25.0, 'degC']
        :param attrs: Dictionary of attributes and metadata. Entries with 'name', 'long_name', 'title',
            'summary' or others, especially as per netCDF attribute conventions.
            The 'manufacturer' should possibly also be provided.
        :return:
        """
        # Deal with the pitch of the pixels (centre-to-center spacing)
        # TODO : Reconsider storage of scalar quantities with units and attributes
        # TODO : Multiple channels - implies multi-axis SQE/ASR
        # TODO : SQE channel axis ('R', 'G', 'B' ?)
        # TODO : Add name field to FocalPlaneArray as well as long_name
        if len(pitch) == 3:
            self.pitchx = check_convert_units([pitch[0], pitch[2]], 'mm')
            self.pitchy = check_convert_units([pitch[1], pitch[2]], 'mm')
        elif len(pitch) == 2:
            self.pitchx = self.pitchy = check_convert_units(pitch, 'mm')
        else:
            warnings.warn('The input "pitch" to FocalPlaneArray must be 2 or 3 element list with pitchx, '
                          'pitchy (if different) and units.')
        self.pitch_units = 'mm'
        # Deal with pixel aperture
        if len(aperture) == 3:
            self.aperturex = check_convert_units([aperture[0], aperture[2]], 'mm')
            self.aperturey = check_convert_units([aperture[1], aperture[2]], 'mm')
        elif len(aperture) == 2:
            self.aperturex = self.aperturey = check_convert_units(aperture, 'mm')
        else:
            warnings.warn('The input "aperture" to FocalPlaneArray must be 2 or 3 element list with aperturex, '
                          'aperturey (if different) and units.')
        self.pixels = pixels  # number of pixels in x and y directions
        self.wellcapacity = wellcapacity  # should be a scalar value in electrons
        self.readnoise = readnoise  # scalar in electrons
        # Darkcurrent a little more complicated to deal with - need to check if given per unit area or per pixel
        q_darkcurrent = Q_(*darkcurrent)
        if '[length]' in q_darkcurrent.dimensionality.keys():  # The dark current is presumably given per unit area
            q_darkcurrent = q_darkcurrent.to('e/s/mm^2')
            q_darkcurrent *= Q_(self.pitchx * self.pitchy, 'mm^2')
            darkcurrent = [q_darkcurrent.magnitude, 'e/s']
        self.darkcurrent = check_convert_units(darkcurrent, 'e/s')  # Actually e/s/pixel
        self.darkcurrent_ref = self.darkcurrent
        self.darkcurrent_units = 'e/s'
        self.dsnu = dsnu
        self.prnu = prnu
        # Deal with ASR or SQE
        if (sqe is None) and (asr is None):
            warnings.warn('Either SQE or ASR (but not both) must be provided for FocalPlaneArray objects.')
        elif (sqe is not None) and (asr is not None):
            warnings.warn('Either SQE or ASR (but not both) must be provided for FocalPlaneArray objects.')
        elif asr is not None:  # Convert to SQE as primary required quantity
            xd_check_convert_units(asr, 'wvl', 'nm')
            xd_check_convert_units(asr, 'asr', 'A/W')
            self.asr = asr
            self.sqe = xd_asr2sqe(asr)
        else:
            xd_check_convert_units(sqe, 'wvl', 'nm')
            self.sqe = sqe
            self.asr = xd_sqe2asr(sqe)
        self.t_ref = check_convert_units(t_ref, 'degC')
        self.temperature_units = 'degC'
        self.darkcurrent_delta_t = check_convert_units(darkcurrent_delta_t, 'delta_degC')
        self.temperature = temperature  # Will also set dark current
        self.attrs = attrs  # Attach user-defined attributes
        # Calculate the horizontal and vertical MTF of the array
        # First set up a set of spatial frequencies up to a factor of the pixel nyquist
        nyquist_x = 1.0/(2.0*self.pitchx)  # cy/mm
        self.nyquist_x = nyquist_x
        nyquist_y = 1.0/(2.0*self.pitchy)  # cy/mm
        self.nyquist_y = nyquist_y
        nyquist_max = np.maximum(nyquist_x, nyquist_y)
        nyquist_min = np.minimum(nyquist_x, nyquist_y)
        # Generate a relative set of spatial frequencies, with variable spacing up to 20 times relative nyquist
        spf_rel = np.hstack((np.linspace(0., 1, 10), np.linspace(1., 2, 9), np.linspace(2. ,3, 7),
                             np.linspace(3., 4, 7),  np.linspace(4., 5, 7), np.linspace(5. ,6, 5),
                             np.linspace(6., 7, 5),  np.linspace(7., 8, 5), np.linspace(8. ,9 ,3),
                             np.linspace(9., 10, 3)))  # Will use sinc function up to argument of 10
        spf_rel = np.unique(spf_rel)   # Will be duplication at the zeroes
        spf_x = spf_rel * nyquist_x * 2.0
        spf_y = spf_rel * nyquist_y * 2.0
        spf = xd_identity(np.unique(np.hstack((spf_x, spf_y))), 'spf', attrs={'units': '1/mm'})  # Create spat freq axis
        self.spf = spf
        # Get back spf_rel, which could be different in x and y directions
        spf_rel_x = spf.data / (nyquist_x * 2.0)
        spf_rel_y = spf.data / (nyquist_y * 2.0)
        fldo = xd_identity([0.0, 90.0], 'fldo')
        self.mtf = xray.DataArray(np.sinc(np.vstack((spf_rel_x, spf_rel_y))).T,
                                  [(spf), (fldo)],
                                  name='mtf', attrs={'units': ''})


    def set_dark_current(self):
        """ Set the dark current of the FocalPlaneArray (FPA) according to the dark current reference temperature and
            the current operating temperature of the FPA.
        :return:
        """
        temp_difference = self.temperature - self.t_ref
        factor = 2.0**(temp_difference / self.darkcurrent_delta_t)
        self.darkcurrent = self.darkcurrent_ref * factor
