__author__ = 'DGriffith, ARamkolowan'

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
from .. import ureg, Q_, U_  # Import the pint units registry from parent
from ..tools.xd import *  # Import additional tools for working with xray.DataArray objects
from ..moglo import *  # Import global glossary/vocalbulary

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


class FocalPlaneArray():
    """ Focal plane array detector. This implementation is typically at the chip level. That is, all or most of the
        information for building an FPA object can be found in the chip-level datasheet. The FPAs in question here
        are usually CCD, CMOS, scientific CMOS or electron-multiplying CCD (EMCCD).

        The FPA class can be combined with an image intensifier tube (IIT) to produce an ICCD device.

        This class does not model Time-Delay and Integration (TDI), which is a dynamic imaging process.


    """

    def __init__(self, pitch, aperture, wellcapacity, readnoise, darkcurrent, dsnu, prnu, sqe=None, asr=None):
        """ Constructor for FocalPlaneAArray objects

        :param pitch: The centre-to-centre spacing of the FPA detector elements. This must be a list where the
            first element is the centre-to-centre spacing in the x-direction (or both x and y), the second element is
            the spacing in the y-direction (if different from the x-direction) and the last element is the units
            in which the pitch is provided (string). The units must be pint-recognizable.
        :param aperture: The effective pixel aperture of the FPA detector elements. This must also be a list, providing
            the x-aperture of the pixel, y-aperture (if different from x) and the units (string).
        :param wellcapacity: The well capacity of the pixel in number of electrons
        :param readnoise: The RMS read noise in electrons.
        :param darkcurrent: The dark current as a list, giving nagnitude and units. Units can be electrons per pixel
            per second (e/s), or in an amperage per pixel (A/s) or an amperage per unit area of the FPA (e.g. A/cm^2)
        :param dsnu: The dark signal non uniformity provided as a list giving the magnitude in the first element
            and units in the second element. The units can also be e/s, A/s or A/area. This is the standard deviation.
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
        :return:
        """
        # Deal with the pitch of the pixels (centre-to-center spacing)
        if len(pitch) == 3:
            self.pitchx = check_convert_units([pitch[0], pitch[2]], 'mm')
            self.pitchy = check_convert_units([pitch[1], pitch[2]], 'mm')
        elif len(pitch) == 2:
            self.pitchx = self.pitchy = check_convert_units(pitch, 'mm')
        else:
            warnings.warn('The input "pitch" to FocalPlaneArray must be 2 or 3 element list with pitchx, '
                          'pitchy (if different) and units.')
        # Deal with pixel aperture
        if len(aperture) == 3:
            self.aperturex = check_convert_units([aperture[0], aperture[2]], 'mm')
            self.aperturey = check_convert_units([aperture[1], aperture[2]], 'mm')
        elif len(aperture) == 2:
            self.aperturex = self.aperturey = check_convert_units(aperture, 'mm')
        else:
            warnings.warn('The input "aperture" to FocalPlaneArray must be 2 or 3 element list with aperturex, '
                          'aperturey (if different) and units.')
        self.wellcapacity = wellcapacity  # should be a scalar value in electrons
        self.readnoise = readnoise  # scalar in electrons
        self.darkcurrent = check_convert_units(darkcurrent, 'e/s')
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


        # Calculate the horizontal and vertical MTF of the array
