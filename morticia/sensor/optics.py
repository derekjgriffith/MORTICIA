__author__ = 'DGriffith'
__project__ = 'MORTICIA'

"""
.. module:: optics
    :platform: Windows, Unix
    :synopsis:
        The optics module includes all code related to imaging optics as spatial and spectral filters. It also
        includes everything related to light propagation within such imaging optics. It does not include the
        atmospheric radiative transfer code. Functions related to the optical characteristics of the human
        eye are included in this module.

        Some important conventions:
        The most common variable names are:
        wvl : wavelengths, conventionally in nm, but not always, so always make sure
        spf : spatial frequencies, conventionally (for lenses) in cy/mm at the image plane
        fno : focal ratios (ratio of effective focal length to aperture diameter
        wvn : wavenumbers, conventionally in cm^-1

        As a basic check when dealing with wavelengths, the following can be observed
        If the wavelength is:
        < 0.15 : Issue a warning
        > 0.15 and < 150.0 : Assume the spectral variable is wavelength in microns
        > 150.0 and < 15000.0 : Assume the spectral variable is wavelength in nm
        > 15000.0 : Issue a warning

Dependencies : numpy (as np), pandas (as pd) and xray
               easygui, pint, warnings
"""


import numpy as np
import pandas as pd
import xray
import warnings
import logging  # TODO set up global logging if required.


# Import units registry from parent to avoid duplicates
from morticia import ureg, Q_, U_

# Import tools related to xray DataArray and related unit checking/conversion
from morticia.tools.xd import *

# As a general rule, where relevant and present, optical parameters are passed in in the order
#   spatial frequency (spf), wavelength (wvl), focal ratio (fno) followed by any other parameters.
# Returned functions of these parameters are generally multi-dimensional numpy arrays with spatial
# frequencies varying down columns, wavelength varying across rows and focal ratio varying in the depth
# dimension. However, singleton dimensions are squeezed out using the numpy.squeeze() function. For example,
# if spf and fno are both vectors, but wvl is a scalar, then spf will vary down columns as before, but fno will
# vary along rows rather than in the third dimension.

# TODO : Review the above statements. In general it may not be a good idea to remove singleton dimensions, other
# TODO : than trailing singleton dimensions. The real fix is to use xray.DataArray


# Possible strategy for dealing with units.
# Import pint and define function for creating quantity numpy arrays (typically Q_)
# Alternatively import a "stripped down", locally written version, which provides access to units, but no automatic
# unit checking or conversion


def mtf(spf, wvl, fno):
    """ Compute the simple (optimally focussed) diffraction Modulation Transfer Function (MTF) of a prefect lens with an
    unobscured circular aperture.

    :param spf: Spatial frequencies in the image at which to compute the MTF
    :param wvl: Wavelength in units consistent with the spatial frequencies f
    :param fno: Focal ratio (working focal ratio) of the lens
    :return: Modulation Transfer Function, with spatial frequency (spf) varying down columns and wavelength across rows
        If the frequencies are given in cycles per millimetre, the wavelengths must be in mm.

    Any of the inputs can be a vector. The spatial frequencies are assigned to the rows of the output array, the
    wavelengths vary from column to column and Fno will vary in the third dimension, but singleton dimensions will
    be squeezed out.

    .. seealso:: optics.pmtf, optics.pmtf_obs, optics.pmtf_obs_wfe
    """
    wvl, spf, fno = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(), np.asarray(spf, dtype=np.float64).ravel(),
                                np.asarray(fno, dtype=np.float64).ravel())
    # Compute the cutoff frequencies
    cutoff = 1.0 / (wvl * fno)
    # Any spatial frequencies above the cutoff are set to the cutoff frequency
    spf = np.minimum(spf, cutoff)
    phi = np.arccos(fno * spf * wvl)
    csphi = np.cos(phi) * np.sin(phi)
    the_mtf = 2.0 * (phi - csphi) / np.pi
    return the_mtf.squeeze()


def autocorr_circle(circle_radius, shift):
    """ Compute autocorrelation of a circular aperture of radius **circle_radius** with centre-to-centre sample
    displacements of **shift**.

    :param circle_radius: Radius of circle
    :param shift: Centre-to-centre displacements at which to compute the autocorrelation
    :return: Autocorrelation magnitude

    .. seealso:: optics.crosscorr_circle
    """
    circle_radius = np.asarray(circle_radius, dtype=np.complex128)
    shift = np.asarray(shift, dtype=np.complex128)
    shift = np.minimum(shift, 2.0*circle_radius)  # Calculation only valid up to shift = 2 * circle_radius
    surd = np.sqrt(shift**2.0 - 4.0 * circle_radius**2.0)
    auto = 2.0 * circle_radius**2.0 * np.log((surd + shift)/(2.0*circle_radius)) - shift*surd/2.0
    return np.abs(auto)


def crosscorr_circle(circle_radius, shift):
    """ Cross correlation of unit circle with circle of radius **circle_radius**.
    Compute the cross-correlation of a circular aperture of unit radius at the origin,
    with a circular aperture of radius **circle_radius**, with centre-to-centre displacements of **shift**.

    :param circle_radius: Radius of circle to cross-correlate with unit circle (scalar numeric)
    :param shift: Centre-to-centre displacements at which to compute the cross-correlation
    :return: Cross-correlation magnitude

    .. seealso:: optics.autocorr_circle
    """
    circle_radius = np.asarray(circle_radius, dtype=np.complex128)
    shift = np.asarray(shift, dtype=np.complex128)
    # Formula only valid for shift between 1-circle_radius and 1+circle_radius
    shift = np.minimum(shift, 1 + circle_radius)
    shift = np.maximum(shift, 1 - circle_radius)
    sigma = np.sqrt((circle_radius + shift + 1.0)*(circle_radius - shift + 1.0)*(circle_radius + shift - 1.0)*(
                        circle_radius - shift - 1.0))
    cross = circle_radius**2.0 * np.log(2.0*circle_radius*shift / (shift**2.0 + circle_radius**2.0 - sigma - 1)) + (
                     np.log((shift**2.0 - circle_radius**2.0 + sigma + 1.0)/(2.0*shift)) - sigma/2.0)
    return np.abs(cross)


def mtf_obs(spf, wvl, fno, obs=0.0):
    """ Compute the optimally focussed diffraction Modulation Transfer Function of a perfect lens with an
    circular aperture having a centred circular obscuration. This is the monochromatic MTF computed at
    a number of discrete given wavelengths.

    :param spf: Spatial frequencies in the image at which to compute the MTF (numpy vector).
    :param wvl: Wavelength in units consistent with the spatial frequencies **spf** (numpy vector)
    :param fno: Focal ratio (working focal ratio) of the lens (numpy vector).
    :param obs: The obscuration ratio (ratio of obscuration diameter to total aperture diameter).
                The obs input must be a scalar numeric.
    :return: MTF with respect to spatial frequency, wavelength, focal ratio and obscuration ratio.
             Singleton dimensions are squeezed out.

    """
    if obs == 0.0:  # just calculate the unobscured diffraction MTF
        return mtf(spf, wvl, fno)

    # Otherwise mesh things up and do it the hard way
    wvl, spf, fno = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(),
                                np.asarray(spf, dtype=np.float64).ravel(),
                                np.asarray(fno, dtype=np.float64).ravel())
    # Calculate w at each matrix site
    w = 2.0 * fno * spf * wvl
    the_mtf = (autocorr_circle(1.0, w) - 2.0*crosscorr_circle(obs, w) + autocorr_circle(obs, w)) / (
                        np.pi*(1.0 - obs**2))
    return np.maximum(the_mtf.squeeze(), 0.0)


def pmtf(spf, wvl, fno, wvl_weights):
    """ Compute polychromatic MTF of lens at given spatial frequencies (in the image plane) for specified wavelengths
    and wavelength weighting factors.

    :param spf: Spatial frequencies in the image at which to compute the polychromatic MTF (numpy vector).
    :param wvl: Wavelength in units consistent with the spatial frequencies **spf** (numpy vector)
    :param fno: Focal ratio (working focal ratio) of the lens (numpy vector).
    :param wvl_weights: A numpy vector having the same length as the wvl vector, providing the relative weights of each
        of the wavelengths.
    :return: Polychromatic MTF with respect to spatial frequency, wavelength and focal ratio.
        Singleton dimensions are squeezed out of the returned numpy array.
    """
    # size of weights MUST be the same size as the wvl
    wvl = np.asarray(wvl, dtype=np.float64).ravel()
    wvl_weights = np.asarray(wvl_weights, dtype=np.float64).ravel()
    spf = np.asarray(spf, dtype=np.float64).ravel()
    fno = np.asarray(fno, dtype=np.float64).ravel()
    if wvl.size != wvl_weights.size:
        logging.error('Number of wavelength weights must be equal to number of wavelengths in call to optics.pmtf()')
    # Tile the weights up to the same size as the meshgridded arrays
    wvl_weights_1 = np.reshape(wvl_weights, (1, wvl_weights.size, 1))
    wvl_weights_2 = np.tile(wvl_weights_1, (spf.size, 1, fno.size))
    wvl, spf, fno = np.meshgrid(wvl, spf, fno)
    # Compute the cutoff frequencies
    cutoff = 1.0 / (wvl * fno)
    # Any spatial frequencies above the cutoff are set to the cutoff frequency
    spf = np.minimum(spf, cutoff)
    phi = np.arccos(fno * spf * wvl)
    csphi = np.cos(phi) * np.sin(phi)
    the_mono_mtf = 2.0 * (phi - csphi) / np.pi
    the_poly_mtf = np.sum(the_mono_mtf * wvl_weights_2, axis=1) / np.sum(wvl_weights_2, axis=1)
    return the_poly_mtf.squeeze()


def pmtf_obs(spf, wvl, fno, wvl_weights, obs=0.0):
    """ Compute polychromatic MTF of obscured lens at given spatial frequencies (in the image plane) for specified
    wavelengths and wavelength weighting factors. The lens may have a circular obcuration of specific ratio.

    :param spf: Spatial frequencies in the image at which to compute the polychromatic MTF (numpy vector).
    :param wvl: Wavelength in units consistent with the spatial frequencies **spf** (numpy vector)
    :param fno: Focal ratio (working focal ratio) of the lens (numpy vector).
    :param wvl_weights: A numpy vector having the same length as the wvl vector, providing the relative weights of each
        of the wavelengths.
    :param obs: Obscuration ratio
    :return: Polychromatic obscured MTF with respect to spatial frequency, wavelength and focal ratio.
        Singleton dimensions are squeezed out of the returned numpy array.
    """
    wvl = np.asarray(wvl, dtype=np.float64).ravel()
    wvl_weights = np.asarray(wvl_weights, dtype=np.float64).ravel()
    spf = np.asarray(spf, dtype=np.float64).ravel()
    fno = np.asarray(fno, dtype=np.float64).ravel()
    if obs == 0.0:  # Just return the unobscured polychromatic MTF
        return pmtf(spf, wvl, fno, wvl_weights)
        # size of weights MUST be the same size as the wvl
    if wvl.size != wvl_weights.size:
        logging.error('Number of wavelength weights must be equal to number of wavelengths in call to optics.pmtf_obs()')
    # Tile the weights up to the same size as the meshgridded arrays
    wvl_weights_1 = np.reshape(wvl_weights, (1, wvl_weights.size, 1))
    wvl_weights_2 = np.tile(wvl_weights_1, (spf.size, 1, fno.size))
    # Calculate the obscured, monochromatic MTFs
    wvl, spf, fno = np.meshgrid(wvl, spf, fno)
    # Calculate w at each matrix site
    w = 2.0 * fno * spf * wvl
    the_mono_mtf = (autocorr_circle(1.0, w) - 2.0*crosscorr_circle(obs, w) + autocorr_circle(obs, w)) / (np.pi*(1.0 - obs**2))
    the_poly_mtf = np.sum(the_mono_mtf * wvl_weights_2, axis=1) / np.sum(wvl_weights_2, axis=1)
    return the_poly_mtf.squeeze()


def atf(spf, wvl, fno, rms_wavefront_error):
    """ Compute the MTF degradation factor for a lens operating at the given wavelengths and and with the given
    focal ratios, RMS wavefront errors at the the specified spatial frequencies in the image plane.
    ATF stands for Aberration Transfer Function.

    :param spf: Spatial frequencies in the image plane at which to compute the ATF. Spatial frequencies must be in
        reciprocal units to wavelengths i.e. if wavelengths are in mm, spatial frequencies must be in cycles per mm.
    :param wvl: Wavelengths at which to compute the ATF
    :param fno: Focal ratios at which to compute the ATF
    :param rms_wavefront_error: RMS wavefront error magnitudes (in waves) at which to compute the ATF
    :return: A numpy array with the aberration transfer function.

    Reference : Shannon, R.R., Handbook of Optics, Volume 1, 2nd Edition, Chapter 35 - Optical
    Specifications. "This is an approximation, however, and it becomes progressively less accurate as
    the amount of the rms wavefront error exceeds about 0.18 wavelength."

    The formula used for computing the aberration MTF due to RMS wavefront error of :math:`W` at spatial frequency
    :math:`f` is

    .. math::
        M\\!T\\!F_{W}(f)=1-\\left(\\frac{W}{0.18}\\right)^{2}\\left[1-4\\left(\\frac{f}{f_{c}}-\\frac{1}{2}\\right)^{2}\\right]

    where the diffraction cutoff (or "critical") frequency is :math:`f_c`.

    .. seealso:: optics.pmtf_obs_wfe
    """
    rms_wavefront_error = np.abs(rms_wavefront_error)  # Force positive
    if np.max(rms_wavefront_error) > 0.3:
        logging.warning('optics.atf function generally only valid up to RMS wavefront error of 0.3. Called with'
                        'maximum value of %f.', np.max(rms_wavefront_error))
    wvl, spf, fno, rms_wavefront_error = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(),
                                                     np.asarray(spf, dtype=np.float64).ravel(),
                                                     np.asarray(fno, dtype=np.float64).ravel(),
                                                     np.asarray(rms_wavefront_error, dtype=np.float64).ravel())
    # Find the cutoff frequencies
    cutoff = 1.0 / (fno * wvl)
    # Compute the spatial frequencies as a fraction of the cutoff
    nu = spf / cutoff
    # Compute the ATF according to Shannon
    the_atf = 1.0 - ((rms_wavefront_error / 0.18)**2.0) * (1.0 - 4.0 * (nu - 0.5)**2.0)
    # Values above 1.0 are not possible, so set those to 1.0
    the_atf[the_atf > 1.0] = 1.0
    return the_atf.squeeze()


def patf(spf, wvl, fno, rms_wavefront_error, wvl_weights):
    """ Compute the polychromatic aberration transfer function.

    :param spf: Spatial frequencies in the image plane at which to compute the ATF. Spatial frequencies must be in
        reciprocal units to wavelengths i.e. if wavelengths are in mm, spatial frequencies must be in cycles per mm.
    :param wvl: Wavelengths at which to compute the ATF
    :param fno: Focal ratios at which to compute the ATF
    :param rms_wavefront_error: RMS wavefront error magnitudes (in waves) at which to compute the ATF
    :param wvl_weights: A numpy vector having the same length as the wvl vector, providing the relative weights of each
        of the wavelengths.
    :return: Polychromatic Aberration Transfer Function in a numpy array.

    .. seealso:: optics.atf
    """
    wvl = np.asarray(wvl, dtype=np.float64).ravel()
    wvl_weights = np.asarray(wvl_weights, dtype=np.float64).ravel()
    spf = np.asarray(spf, dtype=np.float64).ravel()
    fno = np.asarray(fno, dtype=np.float64).ravel()
    rms_wavefront_error = np.asarray(rms_wavefront_error, dtype=np.float64).ravel()
    # size of weights MUST be the same size as the wvl
    if wvl.size != wvl_weights.size:
        logging.error('Number of wavelength weights must be equal to number of wavelengths in call to optics.patf()')
    if wvl.size != rms_wavefront_error.size:
        logging.error('Number of wavefront error elements must be equal to number of wavelengths in'
                      ' call to optics.patf()')
    # Tile the weights up to the same size as the meshgridded arrays
    wvl_weights_1 = np.reshape(wvl_weights, (1, wvl_weights.size, 1))
    wvl_weights_2 = np.tile(wvl_weights_1, (spf.size, 1, fno.size))
    # Also tile up the wavefront errors to the same size
    rms_wavefront_error = np.abs(rms_wavefront_error)  # Force positive
    rms_wavefront_error_1 = np.reshape(rms_wavefront_error, (1, rms_wavefront_error.size, 1))
    rms_wavefront_error_2 = np.tile(rms_wavefront_error_1, (spf.size, 1, fno.size))

    if np.max(rms_wavefront_error) > 0.3:
        logging.warning('optics.atf function generally only valid up to RMS wavefront error of 0.3. Called with'
                        'maximum value of %f.', np.max(rms_wavefront_error))
    wvl, spf, fno = np.meshgrid(wvl, spf, fno)
    # Find the cutoff frequencies
    cutoff = 1.0 / (fno * wvl)
    # Compute the spatial frequencies as a fraction of the cutoff
    nu = spf / cutoff
    # Compute the ATF according to Shannon
    the_mono_atf = 1.0 - ((rms_wavefront_error_2 / 0.18)**2.0) * (1.0 - 4.0 * (nu - 0.5)**2.0)
    the_poly_atf = np.sum(the_mono_atf * wvl_weights_2, axis=1) / np.sum(wvl_weights_2, axis=1)
    # Values above 1.0 are not possible, so set those to 1.0
    the_poly_atf[the_poly_atf > 1.0] = 1.0
    return the_poly_atf.squeeze()



def pmtf_obs_wfe(spf, wvl, fno, rms_wavefront_error, wvl_weights, obs=0.0):
    """ Compute polychromatic modulation transfer function for lens having circular pupil with centred circular
    obscuration and with aberrations expressed in terms of RMS wavefront error.

    :param spf: Spatial frequencies in the image at which to compute the polychromatic MTF (numpy vector).
    :param wvl: Wavelength in units consistent with the spatial frequencies **spf** (numpy vector)
    :param fno: Focal ratio (working focal ratio) of the lens (numpy vector).
    :param rms_wavefront_error: Numpy vector of RMS wavefront error values
    :param wvl_weights: A numpy vector having the same length as the wvl vector, providing the relative weights of each
        of the wavelengths.
    :param obs: Obscuration ratio (centred circular obscuration in circular pupil), ratio of obscuration diameter to
        aperture diameter.
    :return: Numpy array with polychromatic modulation transfer function.

    .. seealso:: optics.atf, optics.pmtf_obs, optics.pmtf
    """
    # Compute the obscured polychromatic MTF
    poly_mtf_obs = pmtf_obs(spf, wvl, fno, wvl_weights, obs)
    # Compute the polychromatic aberration transfer function
    poly_atf = patf(spf, wvl, fno, rms_wavefront_error, wvl_weights)
    # Compute the return result as an element-wise product
    poly_mtf_obs_wfe = poly_mtf_obs * poly_atf
    return poly_mtf_obs_wfe.squeeze()


def n_air(wvl, temperature, pressure):
    """ Return the refractive index of air computed using the same formula used by ZEMAX
    See the section on Index of Refraction Computation in the Thermal Analysis chapter of the ZEMAX manual.

    :param wvl:  Wavelength(s) in microns. If all values of wvl exceed 100, then wavelengths are assumed to be in nm
    :param temperature: Temperature in Celsius.
    :param pressure: Relative air pressure (atmospheres, with 1 atm = 101 325 Pa).
    :return: This function returns a matrix with wvl varying from row to row, temperature varying from column to column
        and pressure varying in the depth dimension. The returned matrix is subject to np.squeeze() to remove any
        singleton dimensions.

    Reference :
    F. Kohlrausch, Praktische Physik, 1968, Vol 1, page 408
    """
    wvl = np.array(wvl, dtype=np.float)
    if np.all(wvl >= 100.0):
        wvl /= 1000.0  # Convert to microns if all wavelengths are greater than 100
    temperature, wvl, pressure = np.meshgrid(temperature, wvl, pressure)
    # Compute the reference refractive indices across all wavelengths
    n_ref = 1. + (6432.8 + (2949810. * wvl**2) / (146. * wvl**2 - 1) + (25540. * wvl**2) / (41. * wvl**2 - 1.)) * 1e-8
    # Compute the full data set (potentially 3D)
    air_rin = 1. + ((n_ref - 1.) * pressure) / (1. + (temperature - 15.) * 3.4785e-3)
    return air_rin.squeeze()


# Functions related to the human eye, namely contrast transfer function (CTF) and modulation transfer function (MTF)
def ctf_eye(spf, lum, w, num_eyes=2, formula=1):
    """ Compute the contrast transfer function of the human eye.
    By default, uses the condensed version of the Barten CTF.

    :param spf: spatial frequencies in eye-space in cycles per milliradian (scalar or vector numpy array input)
    :param lum: mean luminance of the viewing area in :math:`cd/m^2` (scalar or vector numpy array input)
    :param w: the angular width of the viewing area, or the square root of the angular viewing area in square degrees
        (scalar or vector numpy input).
    :param num_eyes: The number of eyes used for viewing (2 for binocular viewing or 1 for monocular viewing). The
        default is num_eyes=2.
    :param formula: The formula variant used for the computation. Defaults (formula=1) to the simple formula first
        published by Barten in SPIE 2003. Other options are formula=11 and formula=14, which are slight variations.
    :return: The CTF with respect to spf, lum and w (up to a 3D numpy array). Singular dimensions are squeezed out
        using numpy.squeeze().
    """
    spf = np.array(spf, dtype=np.float)
    spf, lum, w = np.meshgrid(spf, lum, w)
    # Convert spatial frequencies to cycles per degree
    u = 1000.0 * np.pi * spf / 180.0
    if formula == 1:
        num = (540.0 * (1.0 + 0.7 / lum)**-0.2)
        denom = (1.0 + 12.0 / (w * (1.0 + u / 3.0)))
        c = 0.06
        b = 0.3 * (1.0 + 100.0 / lum)**0.15
        a = num / denom
        thresh = 1.0 / (a * u * np.exp(-b * u) * np.sqrt(1.0 + c * np.exp(b * u)))
    elif formula == 14:
        num = (540.0 * (1.0 + 0.7 / lum)**-0.2)
        denom = (1.0 + 12.0 / (w * (1.0 + u / 3.0)**2))  # Notice square on 1+u/3 factor
        c = 0.06
        b = 0.3 * (1.0 + 100.0 / lum)**0.15
        a = num / denom
        thresh = 1.0 / (a * u * np.exp(-b * u) * np.sqrt(1.0 + c * np.exp(b * u)))
    elif formula == 11:  # A more complex formula that also only uses the viewing size and mean luminance
        m_opt = np.exp(-0.0016 * u**2 * (1.0 + 100.0 / lum)**0.08)
        num = 5200.0 * m_opt
        denom = np.sqrt((1.0 + 144. / w**2 + 0.64 * u**2)*(63. / lum**0.83 + 1.0 / (1.0 - np.exp(-0.02 * u**2))))
        thresh = denom / num
    else:  # create array of nans
        thresh = np.zeros(shape=spf.shape)
        thresh[:] = np.nan
    thresh = np.squeeze(thresh)
    if num_eyes == 1:
        thresh *= np.sqrt(2.0)
    thresh[thresh > 1.0] = np.nan  # Threshold greater than 1 is meaningless
    return thresh



class Lens(object):
    """ The Lens class encapsulates information and behaviour related to imaging lens systems.
    The chief characteristics of a lens are its spectral through-field, through-focus and
    through-frequency MTF, as well as the spectral transmission.
    In order to transform spatial frequencies in the image plane to angular spatial frequencies
    in object space, the effective focal length of the lens (efl) must also be known.
    The basic lens model is a near diffraction-limited system with a centred circular aperture
    having a centred circular obscuration (which may be absent), where the MTF is constant over
    the entire field of view (FOV). A lens with field-dependent MTF can be constructed by
    providing wavefront error input that varies with field.

    The most basic lens model implemented here, from which more complicated lens models could inherit
    their properties have the following attributes

    :param efl: The effective focal length of the lens in mm
    :param fno: The focal ratio of the lens
    :param trn: The spectral transmission of the lens (zero to unity).
    :param wfe: The wavefront error measured in waves. This can be a scalar, assumed the same for
        all wavelengths, or it can be provided as a function of wavelength and/or field position.
    :param obs: The obscuration ratio, being the ratio of the circular obscuration diameter to the
        full circular aperture aperture diameter
    :param mtf: The MTF of the lens. Either the MTF can be provided as a set of measurements or it can be
        computed from efl, fno, obs and wfe

    The lens MTF is computed as a function of spatial frequency in the image, wavelength, defocus and field position.

    The total RMS wavefront error is computed as

    .. math::
        W=\\sqrt{W_{\\Delta\\!z}^{2}+W_{a}^{2}}=\\sqrt{\\left(\\frac{\\Delta\\!z}{8\\lambda F^{2}}\\right)^{2}+W_{a}^{2}}

    where :math:`\\Delta\\!z` is the defocus expressed in the same units as the wavelength :math:`\\lambda`, :math:`F` is
    the focal ratio and :math:`W_a` is the RMS waverfront error due to aberrations at best focus.

    """


    def __init__(self, efl, fno, trn, wfe=None, obs=None, mtf=None, wvn_step=500.0, wfe_allowed=0.3, attrs=None):
        """ Lens constructor.
        The lens is constructed using the focal length, focal ratio and spectral transmittance, and
        optionally also the obscuration ratio and wavefront error.

        :param efl: The effective focal length of the lens. The effective focal length must be a scalar value
            with units, e.g. [30, 'mm'].
        :param fno: The focal ratio (or f-number) which is th ratio of focal length to aperture diameter. This input
            must be a scalar value.
        :param trn: The spectral transmission function, typically a function of wavelength. The spectral
            transmittance function must be
        :type trn: xray.DataArray
        :param wfe: The wavefront error must be expressed in waves (unitless). It can be a function of
            wavelength (wvl) and field position (flr) , but not focus (fldz). If no input is provided, the wfe
            will default to zero over the same spectral region as the lens transmission function.
        :param obs: The obscuration ratio of the lens (scalar numeric) range 0.0 - 1.0, Default None (0.0)
        :param mtf: A complete pre-computed or measured MTF, with axes of spf, wvl, fldz, flr, flo
        :param wvn_step: Minimum spectral increment in wavenumbers (cm^-1) for MTF calculation. Default 500 cm^-1
        :param wfe_allowed: Maximum allowed wfe when computing lens defocus wfe, Default 0.5 waves at mean wavelength
            of transmission function domain.
        :param attrs: A dictionary of user-defined attributes and metadata. Consider including 'name' and
            'long_name', 'title', 'summary' or other netCDF convention attributes.
        :return:
        """

        # Check some assertions : this is bad practice - assertions are used to trap situations
        # the demonstrate that there is a bug in the code. Rather deal with input checking
        # in other ways
        if not trn.dims == ('wvl',):  # TODO : Throw exception, issue warning, log warning
            warnings.warn('Data axes for transmission of optics.Lens are likely missing/incorrect.')
            raise MissingDataArrayAxis('Transmission data for optics.Lens to be function of wavelength only.')
        # Check units of transmission wavelength scale and convert
        xd_check_convert_units(trn, 'wvl', 'nm')  # Change units on wvl axis to nm in place
        if any(trn['wvl'] < 150.0) or any(trn['wvl'] > 15000.0):
            warnings.warn('Wavelength units for optics.Lens transmission probably not in units of '+ trn.attrs['wvl_units'])
        self.trn = trn
        # Check units of efl and convert
        self.efl = Scalar('efl', *efl)  # convert efl units to default units
        self.fno = Scalar('fno', fno, '')  # fno is unitless
        self.trn = trn  # this should be an xray.DataArray
        if obs:
            if obs < 0.0 or obs > 1.0: warnings.warn('Obscuration ratio for optics.Lens must be from 0.0 to 1.0')
            self.obs = Scalar('obs', obs, '')  # Obscuration is unitless
        else:
            self.obs = None  # no obscuration
        # Need to choose the wavelength grid on which to compute the MTF
        # At this point is is assumed that EFL is in mm and wavelengths are in nm
        wvl_min = np.min(trn['wvl']).values
        wvl_max = np.max(trn['wvl']).values
        wvl_mean = np.mean([wvl_min, wvl_max])
        if wvl_min == wvl_max:
            warnings.warn('Spectral transmission function of optics.Lens should span a wavelength region.')
        # Spectral points will be evenly spaced in wavenumber
        wvn_min = 1.0e7 / wvl_max  # cm^-1
        wvn_max = 1.0e7 / wvl_min  # cm^-1
        nsteps = np.ceil((wvn_max - wvn_min) / wvn_step)
        wvn = np.linspace(wvn_max, wvn_min, nsteps+1)  # cm^-1
        wvl = 1.0e7 / wvn  # nm
        # Convert wvl back to identity DataArray
        wvl = xd_identity(wvl, 'wvl')
        self.wvl = wvl
       # Need to choose the spatial frequency grid on which to computer the MTF
        # This depends on the cutoff (or "critical") frequency mainly.
        # Determine the minimum and maximum cutoff frequencies
        cutoff_max = 1.0 / (1.0e-6 * wvl_min * fno)  # wvl assumed now in nm, freq in cy/mm
        self.cutoff_max = Scalar('spfcut', cutoff_max, '1/mm')
        cutoff_min = 1.0 / (1.0e-6 * wvl_max * fno)  # cy/mm
        spf_step = cutoff_min / 12.0  # Want at least 15 samples up to minimum cutoff
        n_steps_spf = cutoff_max * 1.1 / spf_step  # Want sampling up to 1.1 times maximum spf
        spf = np.linspace(0.0, cutoff_max * 1.1, n_steps_spf + 1)  # cy/mm
        # Create a DataArray of spatial frequencies
        spf = xd_identity(spf, 'spf')
        self.spf = spf
        self.spf_max = Scalar('spf', np.max(spf.data), '1/mm')
        # Need to compute the defocus grid on which to compute the MTF
        # There is no real point in computing MTF using the Shannon formula
        # if the total wavefront deformation is greater than 0.18 waves.
        # Will allow up to 0.5 waves of defocus in steps of 0.02 waves
        # Wavefront deformation can be a function of up to 3 variables:
        # wvl, fld and image azimuth (sagittal/tangential)
        # First find the maximum aberration wavefront deformation
        if wfe is not None:
            wfe_max = np.max(wfe)  # waves, RMS
            wfe_wvl_min = np.min(wfe['wvl']).values
            wfe_wvl_max = np.max(wfe['wvl']).values
            if (wfe_wvl_min != wvl_min) or (wfe_wvl_max != wvl_max):
                warnings.warn('Wavefront error input to optics.Lens does not have same wavelength limits as trn.')
        else:
            wfe_max = 0.0
            # Default to zero
            wfe = xray.DataArray([0.0, 0.0], [('wvl', [wvl_min, wvl_max], {'units': 'nm'})], name='wfe',
                                 attrs={'units': ''})
        if wfe_max >= 0.5:
            warnings.warn('WFE for optics.Lens exceeds 0.5 waves. ATF generally invalid for large WFE.')
        defocus_max = 3.5 * 8.0 * 1.0e-6 * wvl_mean * fno**2.0 * np.sqrt(wfe_allowed**2.0 - wfe_max**2.0)
        # Calculate defocus positions in mm
        z_defocus = np.linspace(0.0, defocus_max, 10)  # 10 defocus positions in mm
        # Turn this into an identity DataArray
        z_defocus = xd_identity(z_defocus, 'fldz')
        # Now calculate the wfe values for the different defocus positions, bearing in mind that the wfe could
        # have been given with up to 3 axes (wavelength, field position and field orientation)
        # Axes of wfe should only include wavelength and field position. Those axes must now
        # be complemented with a defocus axis.
        # Calculate the RMS wavefront error due to defocus
        rms_wfe_defocus  = z_defocus / (3.5 * 8.0 * 1e-6 * wvl * fno**2.0)  # Shannon formula, wvl and fldz axes
        rms_wfe_defocus.name = 'wfe'
        rms_wfe_defocus.attrs = {'units': ''}
        # Combine this in quadrature w, but first need to harmonise rms_wfe_defocus with wfe
        rms_wfe_aberr, rms_wfe_defocus = xd_harmonise_interp([wfe, rms_wfe_defocus])
        self.rms_wfe_aberr = rms_wfe_aberr
        self.rms_wfe_defocus = rms_wfe_defocus
        # Calculate total wfe as quadrature sum of aberration and defocus wfe
        self.rms_wfe_total = np.sqrt(self.rms_wfe_aberr**2.0 + self.rms_wfe_defocus**2.0)
        # Set up attributes of the rms_wfe_total
        self.rms_wfe_total.attrs = {}
        self.rms_wfe_total.attrs.update(rms_wfe_aberr.attrs)
        self.rms_wfe_total.attrs.update(rms_wfe_defocus.attrs)
        # Now compute the multidimensional MTF
        self.mtf_obs_wfe()
        # The MTF is a function of spatial frequency, wavelength and defocus (always), with possible extra dimensions
        # of field position (fldx, fldy) and field orientation (fldo).
        self.attrs = attrs  # attach user-defined attribute dictionary


    def mtf_obs_wfe(self):
        """ Compute the multidimensional MTF of a Lens class object

        .. seealso::
            - optics.atf which shows the method and formula by which the aberration/defocus MTF degradation is
                computed.
            - optics.mtf_obs shows documents the typical calculation formulas for obscured MTF
        """
        if self.obs is not None:  # Compute the obscured MTF
            # Calculate w at each matrix site
            w = 2.0 * self.fno.data * self.spf * 1.0e-6 * self.wvl
            # print w
            the_mtf = (autocorr_circle(1.0, w.data) - 2.0*crosscorr_circle(self.obs.data, w.data) +
                       autocorr_circle(self.obs.data, w.data)) / (np.pi*(1.0 - self.obs.data**2.0))
            the_mtf = np.maximum(the_mtf, 0.0)
            # print the_mtf
            # Recreate a data array
            the_mtf = xray.DataArray(the_mtf, [self.spf, self.wvl], name='mtf')
        else:  # Compute the unobscured MTF
            phi = np.arccos(self.fno.data * self.spf * 1.0e-6 * self.wvl)
            csphi = np.cos(phi) * np.sin(phi)
            the_mtf = 2.0 * (phi - csphi) / np.pi
            the_mtf.data[np.isnan(the_mtf.data)] = 0.0  # Set MTF to zero if nan
        # Compute the aberration/defocus transfer function (ATF)
        # Find the cutoff frequencies
        cutoff = 1.0 / (self.fno.data * 1.0e-6 * self.wvl)  # cy/mm
        # Compute the spatial frequencies as a fraction of the cutoff (relative spf)
        nu = self.spf / cutoff
        # Compute the ATF according to Shannon
        the_atf = 1.0 - (1.0 - 4.0 * (nu - 0.5)**2.0) * ((self.rms_wfe_total / 0.18)**2.0)
        # Values above 1.0 are not possible, so set those to 1.0
        the_atf.data[the_atf.data > 1.0] = 1.0
        # Take the product of the MTF and the ATF
        the_mtf = the_mtf * the_atf
        the_mtf.attrs['units'] = ''  # MTF is a unitless quantity
        # Merge Attributes of contributing axes
        the_mtf['spf'].attrs.update(self.spf.attrs)
        the_mtf['wvl'].attrs.update(self.wvl.attrs)
        the_mtf['fldz'].attrs.update({'units': 'mm', 'long_name': long_name['fldz']})
        the_mtf.name = 'mtf'
        self.mtf = the_mtf

    def __repr__(self):
        return (repr(self.efl) + '\n' + repr(self.fno) + '\n' + 'Attributes :' + '\n' +
                repr([name for name in self.__dict__.keys() if not name.startswith('_')]))







