__author__ = 'Ari Ramkilowan'
__project__ = 'MORTICIA'

# My comment for the Merge Test - Ari
# My comment takes two lines.

"""
.. module:: optics
    :platform: Windows, Unix
    :synopsis: The optics module includes all code related to imaging optics as spatial and spectral filters. It also
               includes everything related to light propagation within such imaging optics. It does not include the
               atmospheric radiative transfer code. Functions related to the optical characteristics of the human
               eye are included in this module.
"""

import numpy as np
import logging
# Import global units registry if it exists
# from . import ureg, Q_


def mtf(spf, wvl, fno):
    """
    mtf : Computes the simple (optimally focussed) diffraction Modulation Transfer Function of a prefect lens with an
    unobscured circular aperture
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


def ac_circle(e, w):
    """
    ac_circle(e,w) - Autocorrelation of a circular aperture of radius e
    Computes the autocorrelation of a circular aperture of radius e with centre-to-centre displacements of w
    :param e: Radius of circle
    :param w: Centre-to-centre displacements at which to compute the autocorrelation
    :return: Autocorrelation magnitude
    .. seealso:: optics.cc_circle
    """
    e = np.asarray(e, dtype=np.complex128)
    w = np.asarray(w, dtype=np.complex128)
    w = np.minimum(w, 2.0*e)  # Calculation only valid up to w = 2 * e
    surd = np.sqrt(w**2.0 - 4.0 * e**2.0)
    auto = 2.0 * e**2.0 * np.log((surd + w)/(2.0*e)) - w*surd/2.0
    return np.abs(auto)


def cc_circle(e, w):
    """
    cc_circle(e,w) - Cross correlation of unit circle with circle of radius e.
    Computes the cross-correlation of a circular aperture of unit radius at the origin,
    with a circular aperture of radius e, with centre-to-centre displacements of w.
    :param e: Radius of circle to cross-correlate with unit circle (scalar numeric)
    :param w: Centre-to-centre displacements at which to compute the cross-correlation
    :return: Cross-correlation magnitude
    .. seealso:: optics.ac_circle
    """
    e = np.asarray(e, dtype=np.complex128)
    w = np.asarray(w, dtype=np.complex128)
    # Formula only valid for w between 1-e and 1+e
    w = np.minimum(w, 1 + e)
    w = np.maximum(w, 1 - e)
    sigma = np.sqrt((e + w + 1.0)*(e - w + 1.0)*(e + w - 1.0)*(e - w - 1.0))
    cross = e**2.0 * np.log(2.0*e*w / (w**2.0 + e**2.0 - sigma - 1)) + (
                     np.log((w**2.0 - e**2.0 + sigma + 1.0)/(2.0*w)) - sigma/2.0)
    return np.abs(cross)


def mtf_obs(spf, wvl, fno, obs=0.0):
    """
    mtf_obs : Computes the optimally focussed diffraction Modulation Transfer Function of a prefect lens with an
    circular aperture having a centred circular obscuration
    :param spf: Spatial frequencies in the image at which to compute the MTF
    :param wvl: Wavelength in units consistent with the spatial frequencies f
    :param fno: Focal ratio (working focal ratio) of the lens
    :param obs: The obscuration ratio (ratio of obscuration diameter to total aperture diameter)
                The obs input must be a scalar numeric
    :return: MTF with respect to spatial frequency, wavelength, focal ratio and obscuration ratio
             Singleton dimensions are squeezed out
    """
    if obs == 0.0:  # just calculate the unobscured diffraction MTF
        return mtf(spf, wvl, fno)

    # Otherwise mesh things up and do it the hard way
    wvl, spf, fno = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(), np.asarray(spf, dtype=np.float64).ravel(),
                                np.asarray(fno, dtype=np.float64).ravel())
    # Calculate w at each matrix site
    w = 2.0 * fno * spf * wvl
    the_mtf = (ac_circle(1.0, w) - 2.0*cc_circle(obs, w) + ac_circle(obs, w)) / (np.pi*(1.0 - obs**2))
    return np.maximum(the_mtf.squeeze(), 0.0)


def atf(wvl, fno, rms_wavefront_error, spf):
    """
    atf : Compute the MTF degradation factor for a lens operating at the given wavelengths and and with the given
    focal ratios, RMS wavefront errors at the the specified spatial frequencies in the image plane.
    ATF stands for Aberration Transfer Function.
    :param wvl: Wavelengths at which to compute the ATF
    :param fno: Focal ratios at which to compute the ATF
    :param rms_wavefront_error: RMS wavefront error magnitudes at which to compute the ATF
    :param spf: Spatial frequencies in the image plane at which to compute the ATF
    :return: A numpy array with

    Reference : Shannon, R.R., Handbook of Optics, Volume 1, 2nd Edition, Chapter 35 - Optical
    Specifications. "This is an approximation, however, and it becomes progressively less accurate as
    the amount of the rms wavefront error exceeds about 0.18 wavelength."

    .. seealso:: pmtf_obs_wfe
    """
    if np.max(rms_wavefront_error) > 0.3:
        logging.warning('optics.atf function generally only valid up to RMS wavefront error of 0.3. Called with'
                        'maximum value of %f.', np.max(rms_wavefront_error))


def n_air(wvl, temperature, pressure):
    """
    n_air(wvl, T, P) : Returns the refractive index of air computed using the same formula used by ZEMAX
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
    return np.squeeze(air_rin)


# Functions related to the human eye, namely contrast transfer function (CTF) and modulation transfer function (MTF)
def ctf_eye(spf, lum, w, num_eyes=2, formula=1):
    """
    ctf_eye : The contrast transfer function of the eye.
    By default, uses the condensed version of the Barten CTF

    :param spf: spatial frequencies in eye-space in cycles per milliradian (scalar or vector numpy input)
    :param lum: mean luminance of the viewing area in $cd/m^2$ (scalar or vector numpy input)
    :param w: the angular width of the viewing area, or the square root of the angular viewing area in square degrees
    (scalar or vector numpy input)
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





__author__ = 'ARamkilowan'
