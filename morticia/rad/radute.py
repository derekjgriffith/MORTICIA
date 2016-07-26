__author__ = 'DGriffith, ARamkilowan'

""" This module provides a variety of utilities related to radiometry, including, but not limited to
- Creating, writing and reading of MODTRAN-format .flt (filter) files commonly used to to represent
    spectral power/density functions such as transmittances, reflectance or spectral responsivity functions (SRF).
    The original purpose was to specify sensor spectral response functions (i.e. more than just a "filter").
-
"""

import numpy as np
import pandas as pd
import xarray as xr
import StringIO
import matplotlib.pyplot as plt
import easygui  # for simple file/open dialogs and such
import re
from morticia.tools.xd import *

""" This module provides much of functionality related radiometry required by MORTICIA.
Included here is functionality for :
1) Interfacing to radiative transfer codes.
2) Dealing with spectral filtering and convolution
3) Creation, reading and writing of MODTRAN-style flt filter/SRF definitions

TO DO : Check that MODTRAN can read the .flt files created using the Flt class here
"""

#_micronsymbol = u"\u00B5"
# _micrometres = _micronsymbol + u"m"

# Micron symbol can be encoded in UTF-8 with
_micronsymbol = u'\xB5'.encode('UTF-8')
_micrometres = _micronsymbol + "m"


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


class Flt(object):
    """ Encapsulates a MODTRAN-style .flt spectral response function/filter function definition file.

    """

    def __init__(self, name, units='nm', filterheaders=[], filters=[], centers=[],
                 fwhms=[], shapes=['gauss'], yedges=[0.001], centerflats=[0.0], peakvals=[1.0], wvmins=[], wvmaxs=[],
                 oobs=[np.nextafter(0.0, 1)]):
        """ Create a filter definition object (MODTRAN flt style)
        Input name is mandatory. All other inputs are optional, but the filterheaders list must have the same number of
        string elements as the number of filters. Also, either the filters are given explicitly in the filters input, or
        a list of filter definitions are provided for use with rad.srfgen().
        If not empty, inputs centers through to oobs must be either scalar or have the same number of list elements as
        the filterheaders list. If scalar, the value will be replicated up to the number of filterheader values.

        :param name: Name of the set of filters. If the filterheaders input to this constructor function
            is empty, an attempt will be made to read the data from a file called name + '.flt'
        :param units: Spectral coordinate units for the filter, 'nm', 'um' or 'cm^-1'
        :param filterheaders: List of strings, one header for each filter/SRF in the set.
        :param filters: A list of numpy arrays. Each list element must comprise a 2-column numpy array with the
            spectral coordinate (wavelength in nm or micron or wavenumber per cm) in the first column and the filter
            magnitude in the second column.
        :param centres: rather than provide filters, the inputs to rad.srfgen can be provided, this is a list of center
            wavelengths in nm
        :param fwhms: list of full width at half maximum in nm
        :param shapes: list of strings providing the shapes of the filters (see rad.srfgen for alternatives)
        :param yedges: list of yedge values (see rad.srfgen)
        :param centerflats: list of centerflat values (see rad.srfgen). Note that giving a centerflat value adds this
            amount ot the fwhm of the filters (broadens the resulting width to centerflat + fwhm)
        :param peakvals: list of peak values of the filters
        :param wvmins: list of minimum wavelength limits in nm
        :param wvmaxs: list of maximum wavelength limits in nm
        :param oobs: list of out-of-band leakage values
        :return: MODTRAN-style filter/SRF object
        """
        self.name = name
        if name[-4:].lower() != '.flt':
            self.filename = name + '.flt'
        else:
            self.filename = name
        if not filterheaders:  # Attempt to read from a file
            self.read(self.filename)
            return
        # Check the units input and set up the file header
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
            checkcall = 'use_' + checkitem + ' = Flt.checkparm("' + checkitem + '", ' + checkitem + ', nfilters)'
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
            wvl, y, wvn, wu = srfgen(center=centers[ifilt], fwhm=fwhms[ifilt], shape=shapes[ifilt],
                                 yedge=yedges[ifilt], wvmin=wvmins[ifilt], wvmax=use_wvmaxs[ifilt],
                                 centerflat=centerflats[ifilt], oob=oobs[ifilt], peakval=peakvals[ifilt])
            use_filters.append(np.vstack((wvl, y, wvn, wu)).T)
        self.filters = use_filters


    @staticmethod  # Some input parameter checking for Flt constructor
    def checkparm(parmname, parm, nfilters):
        """ Input parameter checking for Flt constructor

        :param parmname: Name of parameter fpr checking
        :param parm: Parameter value
        :param nfilters: Number of filters
        :return: Checked parameter
        """
        if len(parm) == 1:
            parm *= nfilters
        elif len(parm) != nfilters:
            raise ValueError('Number of ' + parmname + ' must equal number of filterheaders in rad.flt instantiation')
        return parm

    def read(self, filename, name='Unknown'):
        """ Read a .flt format spectral band filter definitions file (MODTRAN format)

        :param filename:
        :return: object of class Flt, if the file is a well-formatted MODTRAN-style .flt file
        """
        isfloatnum = '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'  # regular expression to match a floating point number
        if filename == '.flt':
            filename = easygui.fileopenbox(msg='Please select a .flt file.', filetypes=["*.flt"])
        self.filename = filename
        self.name = name
        with open(filename, 'rt') as fltfil:
            fileheader = fltfil.readline()
            self.unitsheader = fileheader[0].upper()  # Must be W, N or M for wavenumber per cm, nm or microns
            if self.unitsheader == 'W':
                self.units = 'cm^-1'
            elif self.unitsheader == 'N':
                self.units = 'nm'
            elif self.unitsheader == 'M':
                self.units = _micrometres
            else:
                raise ValueError('File header for ' + filename + ' does not start with N, M or W as required for'
                                                                 ' .flt files.')
            self.name = fileheader[1:].strip()
            self.fileheader = fileheader.strip()
            self.filterheaders = []
            self.filters = []
            ifilter = 0  # count the filters
            nexlin = fltfil.readline()
            while nexlin:
                self.filterheaders.append(nexlin.strip())
                nexlin = fltfil.readline()
                if not nexlin:
                    break  # all done
                slin = nexlin.split()  # split line up into a list of tokens at whitespace
                thisfilterdata = np.array([])
                self.filters.append([])
                # if all tokens match
                while all([re.match(isfloatnum, tok) for tok in slin]):  # all tokens a are numbers-this is a line of data
                    # accumulate the data
                    thelinedata = np.array(slin).astype(np.float)
                    if len(thisfilterdata) != 0:
                        thisfilterdata = np.vstack((thisfilterdata, thelinedata))
                    else:
                        thisfilterdata = thelinedata
                    nexlin = fltfil.readline()
                    if not nexlin:
                        break  # all done
                    slin = nexlin.split()  # split line up into a list of tokens at whitespace
                # process the accumulated data depending on what the units are
                if self.unitsheader == 'W':
                    self.filters[ifilter] = np.vstack((1e7/thisfilterdata[:,1], thisfilterdata[:,1],
                                                       thisfilterdata[:,0], 1e9/thisfilterdata[:,1])).T
                elif self.unitsheader == 'N':
                    self.filters[ifilter] = np.vstack((thisfilterdata[:,0], thisfilterdata[:,1],
                                                       1e7/thisfilterdata[:,0], thisfilterdata[:,0]/1000.0)).T
                elif self.unitsheader == 'M':
                    self.filters[ifilter] = np.vstack((thisfilterdata[:,0]*1000.0, thisfilterdata[:,1],
                                                       1e9/thisfilterdata[:,0], thisfilterdata[:,0])).T
                ifilter += 1
            self.nfilters = ifilter


    def __repr__(self, format='  %f'):
        """ Return representation of a .flt MODTRAN-style spectral filter or SRF. This is the same as the format in
         the .flt file is written.

        :param format: format in which to present the numeric data, default is '  %f'
        :return: File representation of .flt object
        """
        selfrep = self.fileheader + '\n'
        for (ifilt, filter) in enumerate(self.filters):
            selfrep = selfrep + self.filterheaders[ifilt] + '\n'
            # Create a string buffer
            strbuff = StringIO.StringIO()
            # column ordering depends on original unit specification
            if self.unitsheader == 'W':
                np.savetxt(strbuff, self.filters[ifilt][:,[2,1,0]], fmt=format)
            elif self.unitsheader == 'N':
                np.savetxt(strbuff, self.filters[ifilt][:,[0,1,2]], fmt=format)
            elif self.unitsheader == 'M':
                np.savetxt(strbuff, self.filters[ifilt][:,[3,1,2]], fmt=format)
            selfrep = selfrep + strbuff.getvalue()
            strbuff.close()  # delete the buffer
        return selfrep

    def write(self, filename=None, format='  %f'):
        """ Write a MODTRAN-style .flt file for this filter/SRF set

        :param filename: Optional filename without extension. If not given, the filter name will be used with
        :param format: Format specifier as for np.savetext for writing the data, default is '  %f'
            extension .flt
        :return: None
        """
        selfrep = self.__repr__(format=format)
        # Write the string to the file
        if filename:
            if filename[-4:].lower != '.flt':
                filename = filename + '.flt'
        else:
            filename = self.filename
        with open(filename, 'wt') as fltfile:
            fltfile.write(selfrep)

    def plot(self):
        """ Plot a MODTRAN-style set of filter/SRF curves.

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
        if len(self.filterheaders) <= 12:
            plt.legend(self.filterheaders, loc='best')
        plt.grid()
        plt.hold(False)

    def flt_as_xd(self):
        """ Convert an Flt class object to a list of xray DataArray objects

        :return: The set of Flt filters as a list of xray DataArray objects, with a wavelength coordinate
            axis ('wvl', long_name = 'Wavelength'
        """
        xd_flt_list = []
        for ifilt in range(self.nfilters):
            wvl = xd_identity(self.filters[ifilt][:,0], 'wvl', self.units)
            #print wvl
            xd_flt_list.append(xr.DataArray(self.filters[ifilt][:,1], [wvl], name='srf',
                                                attrs={'long_name': long_name['srf'],
                                                       'labels': self.filterheaders[ifilt],
                                                       'units': default_units['srf'],
                                                       'title': self.name},
                                                ))
        return xd_flt_list

    def flt_as_xd_harmonised(self, quantity_name='srf', chn_start_index=0):
        """ Convert the Flt class object into a single, wavelength-harmonised xr.DataArray.

        :param quantity_name: The name of the quantity as defined the long_names variable in moglo.py. Defaults to
            'srf', a Spectral Response Function, but could also be a transmission functions 'trn' or other spectral
            quantity known to moglo.py.
        :param chn_start_index: Use this parameter to select the starting channel number. Defaults to zero.
        :return: The set of Flt filters as a single, wavelength-harmonised xr.DataArray object. The filter
            headers are returned in an attribute called 'labels'. The fileheader of the Flt object is
            returned in an attribute called 'title' (netCDF recommendation)
        """
        xd_flt_list = self.flt_as_xd()  # Create a list of xr.DataArray objects
        # Harmonise the wavelength axes
        xd_flt_list_harmonised = xd_harmonise_interp(xd_flt_list)
        xd_attrs_update(xd_flt_list_harmonised)  # Update the attribute
        chn_indices = range(chn_start_index, chn_start_index + len(xd_flt_list))
        # Compile the list into a single object
        flt_data = np.vstack([xd_flt_list_harmonised[ifilt].data for ifilt in range(len(xd_flt_list))])
        xd_flt_harmonised = xr.DataArray(flt_data.T, [(xd_flt_list_harmonised[0]['wvl']),
                                                      ('chn', chn_indices,
                                                         {'labels': self.filterheaders})],
                                           name=quantity_name,
                                           attrs={'long_name': long_name[quantity_name],
                                                  'units': default_units[quantity_name],
                                                  'title': self.name})
        return xd_flt_harmonised

class SpectralChannel(object):
    """ The SpectralChannel class defines any band-limited spectral distribution function. This could be
    the spectral response functions of a sensor, or the spectral transmittance of an optical filter, the spectral
    radiance, irradiance or any other band-limited spectral quantity. The actual spectral distribution
    function is represented by an xr.DataArray object with a wavelength axis in preferred units of 'nm'.

    SpectralChannels is a list of channels that can be indexed in the usual way, by the global channel index
    """
    _channel_counter = 0  # This is a class global counter, incremented for each channel, so that every
                          # defined channel gets a unique number
    _channel_list = []  # The global list of spectral channels indexed self.ichn
    _channel_groups = []  # Global list of channel group names
    _channel_group_dict = {}  # Dictionary of channels indexed by group name
    def __init__(self, sdf=None, group='', channels='all'):
        """ Create one or a group of spectral distribution functions

        :param sdf: The spectral distribution function as an xr.DataArray object.
        :param group:
        :param channels:
        :return:
        """
        pass






