__author__ = 'DGriffith'

import numpy as np
# Import global units registry if it exists
# from . import ureg, Q_



def MTF(spf, wvl, fno):
    '''
    MTF : Computes the simple (optimally focussed) diffraction Modulation Transfer Function of a prefect lens with an
    unobscured circular aperture
    :param spf: Spatial frequencies in the image at which to compute the MTF
    :param wvl: Wavelength in units consistent with the spatial frequencies f
    :param fno: Focal ratio (working focal ratio) of the lens
    :return: Modulation Transfer Function, with spatial frequency (spf) varying down columns and wavelength across rows

    If the frequencies are given in cycles per millimetre, the wavelengths must be in mm.

    Any of the inputs can be a vector. The spatial frequencies are assigned to the rows of the output array, the
    wavelengths vary from column to column and Fno will vary in the third dimension, but singleton dimensions will
    be squeezed out.

    :See Also: PMTF, PMTF, PMTFobs
    '''
    wvl, spf, fno = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(), np.asarray(spf, dtype=np.float64).ravel(),
                                np.asarray(fno, dtype=np.float64).ravel())
    # Compute the cutoff frequencies
    Cutoff = 1.0 / (wvl * fno)
    # Any spatial frequencies above the cutoff are set to the cutoff frequency
    spf = np.minimum(spf, Cutoff)
    phi = np.arccos(fno * spf * wvl)
    csphi = np.cos(phi) * np.sin(phi)
    theMTF = 2.0 * (phi - csphi) / np.pi
    return theMTF.squeeze()

def ACCircle(e, w):
    '''
    ACCircle(e,w) - Autocorrelation of a circular aperture of radius e
    Computes the autocorrelation of a circular aperture of radius e with centre-to-centre displacements of w
    :param e: Radius of circle
    :param w: Centre-to-centre displacements at which to compute the autocorrelation
    :return: Autocorrelation magnitude
    '''
    e = np.asarray(e, dtype=np.complex128)
    w = np.asarray(w, dtype=np.complex128)
    w = np.minimum(w, 2.0*e) # Calculation only valid up to w = 2 * e
    surd = np.sqrt(w**2.0 - 4.0 * e**2.0)
    auto = 2.0 * e**2.0 * np.log((surd + w)/(2.0*e)) - w*surd/2.0
    return np.abs(auto)


def CCCircle(e, w):
    '''
    CCCircle(e,w) - Cross correlation of unit circle with circle of radius e.
    Computes the cross-correlation of a circular aperture of unit radius at the origin,
    with a circular aperture of radius e, with centre-to-centre displacements of w.
    :param e: Radius of circle to cross-correlate with unit circle (scalar numeric)
    :param w: Centre-to-centre displacements at which to compute the cross-correlation
    :return: Cross-correlation magnitude
    '''
    e = np.asarray(e, dtype=np.complex128)
    w = np.asarray(w, dtype=np.complex128)
    # Formula only valid for w between 1-e and 1+e
    w = np.minimum(w, 1 + e)
    w = np.maximum(w, 1 - e)
    sigma = np.sqrt((e + w + 1.0)*(e - w + 1.0)*(e + w - 1.0)*(e - w - 1.0))
    cross = e**2.0 * np.log(2.0*e*w / (w**2.0 + e**2.0 - sigma - 1)) + (
                     np.log((w**2.0 - e**2.0 + sigma + 1.0)/(2.0*w)) - sigma/2.0)
    return np.abs(cross)

def MTFobs(spf, wvl, fno, obs=0.0):
    '''
    MTF : Computes the optimally focussed diffraction Modulation Transfer Function of a prefect lens with an
    circular aperture having a centred circular obscuration
    :param spf: Spatial frequencies in the image at which to compute the MTF
    :param wvl: Wavelength in units consistent with the spatial frequencies f
    :param fno: Focal ratio (working focal ratio) of the lens
    :param obs: The obscuration ratio (ratio of obscuration diameter to total aperture diameter)
                The obs input must be a scalar numeric
    :return:
    '''
    if obs == 0.0:  # just calculate the unobscured diffraction MTF
        return MTF(spf, wvl, fno)

    # Otherwise mesh things up and do it the hard way
    wvl, spf, fno = np.meshgrid(np.asarray(wvl, dtype=np.float64).ravel(), np.asarray(spf, dtype=np.float64).ravel(),
                                np.asarray(fno, dtype=np.float64).ravel())
    # Calculate w at each matrix site
    w = 2.0 * fno * spf * wvl
    theMTF = (ACCircle(1.0, w) - 2.0*CCCircle(obs, w) + ACCircle(obs, w)) / (np.pi*(1.0 - obs**2))
    return np.maximum(theMTF.squeeze(), 0.0)

def n_air(wvl, T, P):
    '''
    n_air(wvl, T, P) : Returns the refractive index of air computed using the same formula used by ZEMAX
    See the section on Index of Refraction Computation in the Thermal Analysis chapter of the ZEMAX manual.

    :param wvl:  Wavelength(s) in microns. If all values of wvl exceed 100, then wavelengths are assumed to be in nm
    :param T: Temperature in Celsius.
    :param P: Relative air pressure (atmospheres, with 1 atm = 101 325 Pa).
    :return: This function returns a matrix with wvl varying from row to row, temperature varying from column to column
    and pressure varying in the depth dimension. The returned matrix is subject to np.squeeze() to remove any
    singleton dimensions.
    Reference :
    F. Kohlrausch, Praktische Physik, 1968, Vol 1, page 408
    '''
    wvl = np.array(wvl, dtype=np.float)
    if np.all(wvl >= 100.0):
        wvl = wvl / 1000.0
    T, wvl, P = np.meshgrid(T, wvl, P)
    # Compute the reference refractive indices across all wavelengths
    n_ref = 1. + (6432.8 + (2949810. * wvl**2) / (146. * wvl**2 - 1) + (25540. * wvl**2) / (41. * wvl**2 - 1.)) * 1e-8
    # Compute the full data set (potentially 3D)
    n_air = 1. + ((n_ref - 1.) * P) / (1. + (T - 15.) * 3.4785e-3)
    return np.squeeze(n_air)
