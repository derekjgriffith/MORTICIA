__author__ = 'DGriffith, ARamkilowan'
import numpy as np
import StringIO
import matplotlib.pyplot as plt

""" This module provides much of functionality related radiometry required by MORTICIA.
Included here is functionality for :
1) Interfacing to radiative transfer codes.
2) Dealing with spectral filtering and convolution
3) Creation, reading and writing of MODTRAN-style flt filter/SRF definitions
"""

_micronsymbol = u"\u00B5"
_micrometres = _micronsymbol + u"m"


def srfgen(center, fwhm, n=101, shape='gauss', yedge=0.001, wvmin=None, wvmax=None, centerflat=0.0,
           oob=np.nextafter(0.0, 1), peakval=1.0, units='nm'):
    """ Generate a spectral filter or spectral response function of various shapes
    :param center: center wavelength of the filter in nm
    :param fwhm: full width at half maximum of the filter  in nm
    :param n: number of spectral samples in the filter (minimum of 3, default 101), should be odd
    :param shape: filter shape, one of 'gauss', 'bartlett', 'welch', 'cosine', 'box', 'cos^2' (default 'gauss')
    :param yedge: minimum y-value at the limits of the filter (default 0.001). No filter values below this threshold
    :param wvmin: extend spectral definition by adding a single point at (wvmin, oob)
    :param wvmax: exptend spectral definition by adding a single point at (wvmax, oob)
    :param centerflat: opens a flat region in the centre of filter having a width of centerflat nm
    :param oob: out-of-band leakage, default 0.0, must be <= yedge
    :param peakval: simply scales the peak of the filter function to this value (default 1.0)
    :param units: spectral axis units, either 'nm', 'cm^-1' or 'um' (nanometers, wavenumber per cm or microns)
    will be returned
    :return: wvlnm, y, wvn, wvlum (wavelengths in nm, filter values, wavenumbers per cm and wavelengths in microns)
    """
    if n < 3:
        raise ValueError('Input parameter n to rad.srfgen must be greater than 3.')
    if oob > yedge:
        raise ValueError('Input parameter oob to rad.srfgen must be smaller than paramter yedge (default 0.001).')
    # force an odd number of points
    if not np.mod(n, 2.0):
        n += 1
    n2 = (n-1.0)/2.0  # mid-point
    x = np.arange(-n2, n2+1, 1.0)  # relative x-coordinates for filter
    shape = shape.lower()
    if shape == 'gauss' or shape == 'gaussian':
        sigma = np.sqrt(-n2**2.0 / 2.0 / np.log(yedge))  # scalar standard deviation
        nfwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma  # normalized fwhm (scalar)
        y = np.exp(-x**2.0 / 2.0 / sigma**2)
    elif shape == 'bartlett' or shape == 'triangle':
        fudge = n2 * (1+yedge+(yedge/100.0))
        nfwhm = n2 + (yedge * fudge)  # normalized fwhm
        y = 1.0 - (abs(x) / nfwhm)  # bartlett
        y[y<0] = np.nextafter(0.0, 1)  # in case of negative samples, make very small positive value
    elif shape == 'welch':
        fudge = n2 * (0.5+(yedge*0.378))
        alph = n2 + (yedge * fudge)
        nfwhm = np.sqrt(2.0) * alph  # normalized fwhm
        y = 1 - ((x ** 2.0) / (alph ** 2.0))  # welch
        y[y<0] = np.nextafter(0.0, 1)  # in case of negative samples, make very small positive value
    elif shape == 'cosine':
        fudge = n2 * (.6365+(yedge*.5))
        alph = n2 + (yedge * fudge)
        nfwhm = (4.0/3.0) * alph  # normalized fwhm
        y = np.cos((np.pi * x) / (2.0 * alph))
        y[y<0] = np.nextafter(0.0, 1)  # in case of negative samples, make very small positive value
    elif shape == 'cos^2' or shape == 'cosquared':
        fudge = n2 * (.6365+(yedge*.5))
        alph = n2 + (yedge * fudge)  # calculate alpha
        nfwhm = alph  # normalized fwhm
        y = np.cos((np.pi * x) / (2.0 * alph))**2
    elif shape == 'box' or shape == 'tophat':
        y = np.ones(x.shape)
        y[0] = np.nextafter(0.0, 1)
        y[-1] = np.nextafter(0.0, 1)
        nfwhm = np.trapz(y)
    else:
        raise ValueError('Unknown SRF/filter shape input to rad.srfgen')
    delta = fwhm / nfwhm  # Determine sample delta (scalar)
    y[y<yedge] = oob  # suppress values dropping below the edge threshold
    wvl = (np.ones(x.shape, dtype=np.float) * center) + (delta * x)
    if centerflat:  # Open a central flat region simply by shifting the wvl-coordinates
        wvl = np.hstack((wvl[:n2+1]-centerflat/2.0, wvl[n2:]+centerflat/2.0))
        y = np.hstack((y[:n2+1], y[n2:]))
    if wvmin and wvmin < wvl[0]:  # insert a point at lower wavelength extremity
        wvl = np.hstack((wvmin, wvl))
        y = np.hstack((oob, y))
    if wvmax and wvmax > wvl[-1]:  # insert a point at upper wavelength extremity
        wvl = np.hstack((wvl, wvmax))
        y = np.hstack((y, oob))
    if units == 'cm^-1':  # wavenumber per cm
        wvn = wvl
        wvlnm = 1.0e7 / wvn
        wvlum = wvlnm / 1000.0
    elif units == 'nm':  # wavelength in nm
        wvlnm = wvl
        wvn = 1.0e7 / wvlnm
        wvlum = wvlnm / 1000.0
    elif units == 'um' or units == _micrometres:  # wavelength in microns
        wvlum = wvl
        wvlnm = 1000.0 * wvlum
        wvn = 1e7 / wvlnm
    else:
        raise ValueError('Wavelength/wavenumber units for MODTRAN .flt definitions must be "cm^1", "nm" or "um".')
    y *= peakval  # finally, scale the peak value
    return wvlnm, y, wvn, wvlum


def tophat(center, fwhm, delta=0.0, wvmin=None, wvmax=None, oob=0.0, units='nm'):
    """ Return tophat/box filter defined by 6 points
    Can also specify out-of-band values and extreme limits, which adds upper and lower bound points
    :param center: center wavelength in nm
    :param fwhm: full width at half max in nm
    :param delta: smallest x-coordinate increment, default see np.nextafter
    :param wvmin: minimum wavelength to reach default None
    :param wvmax: maximum wavelength to reach default None
    :param oob: out-of-band leakage
    :param units: spectral axis units, either 'nm', 'cm^-1' or 'um' (nanometers, wavenumber per cm or microns)
    will be returned
    :return: wvlnm, y, wvn, wvlum (wavelengths in nm, filter values, wavenumbers per cm and wavelengths in microns)
    """
    edgelo = center - fwhm/2.0
    edgehi = center + fwhm/2.0
    if delta != 0.0:
        wvl = np.array([edgelo-delta, edgelo, edgelo+delta, edgehi-delta, edgehi, edgehi+delta])
    else:
        wvl = np.array([np.nextafter(edgelo, -1), edgelo, np.nextafter(edgelo, 1),
                        np.nextafter(edgehi, -1), edgehi, np.nextafter(edgehi, 1)])
    y = np.array([oob, 0.5, 1.0, 1.0, 0.5, oob])
    if wvmin:
        wvl = np.hstack((wvmin, wvl))
        y = np.hstack((oob, y))
    if wvmax:
        wvl = np.hstack((wvl, wvmax))
        y = np.hstack((y, oob))
    if units == 'cm^-1':  # wavenumber per cm
        wvn = wvl
        wvlnm = 1.0e7 / wvn
        wvlum = wvlnm / 1000.0
    elif units == 'nm':  # wavelength in nm
        wvlnm = wvl
        wvn = 1.0e7 / wvlnm
        wvlum = wvlnm / 1000.0
    elif units == 'um' or units == _micrometres:  # wavelength in microns
        wvlum = wvl
        wvlnm = 1000.0 * wvlum
        wvn = 1e7 / wvlnm
    else:
        raise ValueError('Wavelength/wavenumber units for MODTRAN .flt definitions must be "cm^1", "nm" or "um".')
    return wvl, y, wvn, wvlum


class flt:
    @staticmethod
    def checkparm(parmname, parm, nfilters):
        if len(parm) == 1:
            parm *= nfilters
        elif len(parm) != nfilters:
            raise ValueError('Number of ' + parmname + ' must equal number of filterheaders in rad.flt instantiation')
        return parm

    def __init__(self, name, units='nm', filterheaders=[], filters=[], centers=[],
                 fwhms=[], shapes=['gauss'], yedges=[0.001], centerflats=[0.0], peakvals=[1.0], wvmins=[], wvmaxs=[],
                 oobs=[np.nextafter(0.0, 1)]):
        """ Create a filter definition object (MODTRAN flt style)
        Input name is mandatory. All other inputs are optional, but the filterheaders list must have the same number of
        string elements as the number of filters. Also, either the filters are given explicitly in the filters input, or
        a list of filter definitions are provided for use with rad.srfgen().
        If not empty, inputs centers through to oobs must be either scalar or have the same number of list elements as
        the filterheaders list. If scalar, the value will be replicated up to the number of filterheader values.
        :param name: Name of the set of filters
        :param units: Spectral coordinate units for the filter, 'nm', 'um' or 'cm^-1'
        :param filterheaders: List of strings, one header for each filter/SRF in the set.
        :param filters: A list of numpy arrays. Each list element must comprise a 2-column numpy array with the
        spectral coordinate (wavelength in nm ire micron or wavenumer per cm) in the first column and the filter
        magnitude in the second column.
        :param centres: rather than provide filters, the inputs to rad.srfgen can be provided, this is a list of center
        wavelengths in nm
        :param fwhms: list of full width at half maximum in nm
        :param shapes: list of strings providing the shapes of the filters (see rad.srfgen for alternatives)
        :param yedges: list of yedge values (see rad.srfgen)
        :param centerflats: list of centerflat values (see rad.srfgen)
        :param peakvals: list of peak values of the filters
        :param wvmins: list of minimum wavelength limits in nm
        :param wvmaxs: list of maximum wavelength limits in nm
        :param oobs: list of out-of-band leakage values
        :param

        :return: MODTRAN-style filter/SRF object
        """
        self.name = name
        self.filename = name + '.flt'
        if units == 'cm^-1':
            self.unitsheader = 'W' # wavenumber
        elif units == 'nm':
            self.unitsheader = 'N'
        elif units == 'um' or units == _micrometres:
            self.unitsheader = 'M'
        else:
            raise ValueError('Wavelength/wavenumber units for MODTRAN .flt definitions must be "cm^1", "nm" or "um".')
        self.units = units
        self.fileheader = self.unitsheader + ' ' + self.name
        #if len(filterheaders) != len(filters):
        #    raise ValueError('Number of filterheaders must equal number of filters when creating .flt objects.')
        self.filterheaders = filterheaders
        if filters and centers:
            raise ValueError('Do not give both explicit filter definitions and rad.srfgen definitions for rad.flt.')
        nfilters = len(filterheaders)
        self.nfilters = nfilters
        if filters:
            if len(filters) != nfilters:
                raise ValueError('Number of filters must equal number of filter headers in rad.flt instantiation.')
            self.filters = filters  # now fix the filters into the correct column order
            for ifilt in range(nfilters):
                if self.unitsheader == 'W':  # wavenumber
                    # column order is wvnm, y, wvn, wvum
                    self.filters[ifilt] = np.vstack((1.0e7/filters[ifilt][:,0], filters[ifilt][:,1],
                                                     filters[ifilt][:,0], 1.0e7/filters[ifilt][:,0] / 1000.0)).T
                elif self.unitsheader == 'N':
                    self.filters[ifilt] = np.vstack((filters[ifilt][:,0], filters[ifilt][:,1],
                                                     1.0e7 / filters[ifilt][:,0], filters[ifilt][:,0]/1000.0)).T
                elif self.unitsheader == 'M':
                    self.filters[ifilt] = np.vstack((1000.0*filters[ifilt][:,0], filters[ifilt][:,1],
                                                     1.0e10 /filters[ifilt][:,0], filters[ifilt][:,0])).T
            return  # filter definitions have been given explicitly
        # Check that parameters are scalar or multiply them up to size
        checklist = ['centers', 'fwhms', 'shapes', 'yedges', 'centerflats', 'peakvals', 'oobs']
        for checkitem in checklist:
            checkcall = 'use_' + checkitem + ' = flt.checkparm("' + checkitem + '", ' + checkitem + ', nfilters)'
            exec(checkcall)
        if wvmins and len(wvmins) == 1:
            use_wvmins = wvmins * nfilters
        elif wvmins and len(wvmins) != nfilters:
            raise ValueError('Number of wvmins must equal number of filterheaders in rad.flt instantiation')
        else:
            use_wvmins = [[] for ie in range(nfilters)]
        # Check wvmaxs
        use_filters = []
        if wvmaxs and len(wvmaxs) == 1:
            use_wvmaxs = wvmaxs * nfilters
        elif wvmaxs and len(wvmaxs) != nfilters:
            raise ValueError('Number of wvmaxs must equal number of filterheaders in rad.flt instantiation')
        else:
            use_wvmaxs = [[] for ie in range(nfilters)]
        for ifilt in range(nfilters):
            wvl, y, wvn, wu = srfgen(center=use_centers[ifilt], fwhm=use_fwhms[ifilt], shape=use_shapes[ifilt],
                                 yedge=use_yedges[ifilt], wvmin=use_wvmins[ifilt], wvmax=use_wvmaxs[ifilt],
                                 centerflat=use_centerflats[ifilt], oob=use_oobs[ifilt], peakval=use_peakvals[ifilt])
            use_filters.append(np.vstack((wvl, y, wvn, wu)).T)
        self.filters = use_filters

        # Convert the wavelengths to microns or wavenumbers

    def read(filename):
        """ Read a .flt format spectral band filter definitions file (MODTRAN format)
        :param filename:
        :return:
        """
    def __repr__(self):
        """ Return representation of a .flt MODTRAN-style spectral filter or SRF. This is the same as the format in
         the .flt file is written.
        :return: File representation of .flt object
        """
        selfrep = self.fileheader + '\n'
        for (ifilt, filter) in enumerate(self.filters):
            selfrep = selfrep + self.filterheaders[ifilt] + '\n'
            # Create a string buffer
            strbuff = StringIO.StringIO()
            # column ordering depends on original unit specification
            if self.unitsheader == 'W':
                np.savetxt(strbuff, self.filters[ifilt][:,[2,1,0]])
            elif self.unitsheader == 'N':
                np.savetxt(strbuff, self.filters[ifilt][:,[0,1,2]])
            elif self.unitsheader == 'M':
                np.savetxt(strbuff, self.filters[ifilt][:,[3,1,2]])
            selfrep = selfrep + strbuff.getvalue()
            strbuff.close()  # delete the buffer
        return selfrep

    def write(self, filename=None):
        """ Write a MODTRAN-style .flt file for this filter/SRF set
        :param filename: Optional filename without extension. If not given, the filter name will be used with
        extension .flt
        :return: None
        """
        selfrep = self.__repr__()
        # Write the string to the file
        if filename:
            filename = filename + '.flt'
        else:
            filename = self.filename
        with open(filename, 'wt') as fltfile:
            fltfile.write(selfrep)

    def plot(self):
        """ Plot a MODTRAN-style set of filter/SRF curves
        :return: None
        """
        plt.hold(True)
        for ifilt in range(self.nfilters):
            if self.unitsheader == 'W':
                plt.plot(self.filters[ifilt][:, 2], self.filters[ifilt][:,1])
            elif self.unitsheader == 'N':
                plt.plot(self.filters[ifilt][:, 0], self.filters[ifilt][:,1])
            elif self.unitsheader == 'M':
                plt.plot(self.filters[ifilt][:, 3], self.filters[ifilt][:,1])
        plt.title(self.name)
        plt.ylabel('Spectral Response/Transmission')
        if self.units == 'cm^-1':
            plt.xlabel('Wavenumber [' + self.units + ']')
        else:
            plt.xlabel('Wavelength [' + self.units + ']')
        plt.legend(self.filterheaders, loc=0)
        plt.grid()
        plt.hold(False)




