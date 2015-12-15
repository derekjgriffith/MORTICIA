__author__ = 'DGriffith'
"""
 *--------------------------------------------------------------------
 *
 * This file is part of MORTICIA.
 * Copyright (c) 2015-2016 by Derek Griffith and Ari Ramkiloawn
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------"""

"""
.. module:: librad
    :platform: Unix only if running of libRadtran/uvspec is required. Also Windows for setting up uvspec cases.
    :synopsis: This module provides access to libRadtran through the uvspec utility. For information on libRadtran go to
    http://www.libtradtran.org
Some of the code within this module and code imported by the module is provided with libRadtran in the src_py
directory. In order to run libRadtran cases it is necessary to have libRadtran installed on your computer.
There is no Windows-native version of libRadtran, so that generally implies that you are running on a Unix/Linux,
however, reading and writing of uvspec input files, as well as reading of uvspec output files is possible on any
platform.

The relevant code here is taken from libRadtran version 2.0

"""

_isfloatnum = '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'  # regular expression for matching tokens to floating point

import writeLex  # This imports all the libradtran option definitions
import os
import easygui  # For file open dialogs
import numpy as np
import re

uvsOptions = writeLex.loadOptions()  # Load all options into a global dictionary of option specifications.

class Case():
    """ Class which encapsulates a run case of libRadtran/uvspec.
    This class has methods to read libRadtran/uvspec input files, write uvspec input files, run uvspec in parallel on
    multiple compute nodes and read uvspec output files. An important use-case is that of reading a uvspec input
    file called the "base case", altering the parameters of particular option keywords and then running the case
    and reading the outputs. This class is also used by the RadEnv class which encapsulates a radiant environment.
    Construction of radiant environment maps typically requires running an array of librad.Case instances.
    """
    # Definitions of some of the possible uvspec output variables
    uvspecOutVars = {
        'lambda': 'Output wavelengths [nm]',  # Cannot be used in Python because lambda is a keyword use 'wvl
        'wavenumber': 'Wavenumber [cm^1]',  # abbreviate to wvn
        'wvl': 'Short name for wavelengths [nm]',
        'edir': 'Direct beam irradiance (same unit as extraterrestrial irradiance, e.g mW/m^2/nm if using the "atlas3" spectrum in the /data/solar_flux/ directory.)',
        'edn': 'Diffuse downwelling irradiance, i.e. total minus direct beam (same unit as edir).',
        'eup': 'Diffuse upwelling irradiance (same unit as edir).',
        'uavg': 'The mean intensity. Proportional to the actinic flux: To obtain the actinic flux, multiply the mean intensity by 4 pi (same unit as edir).',
        'uavgdir': 'Direct beam contribution to the mean intensity. (same unit as edir).',
        'uavgdn': 'Diffuse downward radiation contribution to the mean intensity. (same unit as edir).',
        'uavgup': 'Diffuse upward radiation contribution to the mean intensity. (same unit as edir).',
        'umu': 'Cosine of the zenith angles for sightline radiance (intensity) calculations.',
        'u0u': 'The azimuthally averaged intensity at numu user specified angles umu. (units of e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)',
        'uu': 'The radiance (intensity) at umu and phi user specified angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)',
        'uu_down': 'The downwelling radiances (intensity) at cmu and phi angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux? directory.)',
        'uu_up': 'The upwelling radiances (intensity) at cmu and phi angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)',
        'cmu': 'Computational polar angles from polradtran',
        'down_fluxI': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes I component). Same units as extraterrestrial irradiance.',
        'up_fluxI': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes I component). Same units as extraterrestrial irradiance.',
        'down_fluxQ': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes Q component). Same units as extraterrestrial irradiance.',
        'up_fluxQ': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes Q component). Same units as extraterrestrial irradiance.',
        'down_fluxU': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes U component). Same units as extraterrestrial irradiance.',
        'up_fluxU': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes U component). Same units as extraterrestrial irradiance',
        'down_fluxV': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes V component). Same units as extraterrestrial irradiance.',
        'up_fluxV': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes V component). Same units as extraterrestrial irradiance.',
        'eglo': 'Global downwelling irradiance',
        'enet': 'Global downwelling minus upwelling (net downward irradiance)',
        'esum': 'Global downwelling plus upwelling',
        'fdir': 'Direct actinic flux (scalar irradiance)',
        'fglo': 'Global actinic flux (scalar irradiance)',
        'fdn': 'Downwelling actinic flux (scalar irradiance)',
        'fup': 'Upwelling actinic flux (scalar irradiance)',
        'uavgglo': 'Total (global) mean diffuse intensity (radiance) = actinic flux/4pi',
        'f': 'Actinic flux (scalar irradiance)',
        'sza': 'Solar zenith angle [deg]',
        'zout': 'Output altitude in km',
        'albedo': 'Albedo',
        'heat': 'Heating rate in K/day'
    }

    # Define a regexp for determining if a token of the zout keyword is a single level
    # It is a single level if it is either a floating point number or the words sur, cpt or toa (for surface,
    # cold point tropopause and top of atmosphere)
    re_isSingleOutputLevel = '(^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$)|(^sur$)|(^cpt$)|(^toa$)'

    def __init__(self, casename='', filename=None, optiondict=None):
        """ A libRadtran/uvspec case
        :param casename: A user-defined name for the libRadtran/uvspec case
        :param filename: An optional filename from which to read the libRadtran/uvspec input
        :param optionlist: A list of option keywords and parameteres (tokens). The keyword existence
         is verified. Besides that, no error checking is performed automatically (yet).
        :return:
        """
        self.name = casename
        self.error_txt = []
        self.options = []  # options is a list [option_name (string), option_tokens (list of strings),
        self.tokens = []
        self.optionobj = []
        self.filorigin = []
        self.solver = 'disort'  # default, modified b y the rte_solver option keyword
        self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
        self.wvl = []  # wavelengths, wavenumbers and output levels difficult to ascertain to start with
        self.wvn = []
        self.zout = []
        self.output_user = ''  # set with the output_user keyword
        self.output_quantity = '' # default is radiances and irradiances
        self.n_umu = 0
        self.umu = []  # zenith angles for radiance calculations
        self.n_phi = 0
        self.phi = []  # azimuth angles for radiance calculations
        self.n_zout = 1  # Assume only one output level, unless zout keyword is used.
        self.n_wvl = '?'  # number of wavelengths is difficult to ascertain to start with
        self.nstokes = 1  # default number of stokes parameters for polradtran
        if filename is not None:
            if not filename:
                # Open a dialog to get the filename
                filename = easygui.fileopenbox(msg='Please select a uvspec input file.', filetypes=["*.INP"])
            elif filename[-4:].lower() != '.inp':
                filename += '.INP'
            opdata, line_nos, path = Case.read(path=filename)
            self.infile = path
            self.outfile = self.infile[:-4] + '.OUT'
            # option_object, source_file_nos (filename, line_number)]
            # Process the results into lists of options
            for (optnumber, option) in enumerate(opdata):
                self.options.append(option[0])  # the option keyword (string)
                self.tokens.append(option[1:])  # The tokens following the keyword (list of strings)
                self.filorigin.append(line_nos[optnumber])  # The file origin of this keyword
                self.optionobj.append(uvsOptions[self.options[optnumber]])  # The option object (dict lookup)
                # Make any possible preparations for occurance of this keyword
                self.prepare_for(option[0],option[1:])
        #TODO Build the case from the dictionary ?


    def prepare_for(self, keyword, tokens):
        """ Make any possible preparations for occurrences of particular keywords

        :param keyword: The uvspec option keyword (string)
        :param tokens: The parameters (tokens) for the keyword as a list of strings
        :return:
        """
        # Prepare for different output formats, depending on the solver
        if keyword == 'rte_solver':
            self.solver = tokens[0]  # First token gives the solver
            if any([self.solver == thesolver for thesolver in ['disort', 'disort2', 'sdisort',
                                                               'spsdisort', 'fdisort1', 'fdisort2']]):
                self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
            elif any([self.solver == thesolver for thesolver in ['twostr','rodents']]):
                self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavg']
            elif self.solver == 'polradtran':
                self.prepare_for_polradtran()
            elif self.solver == 'sslidar':
                self.fluxline = ['center_of_range', 'number_of_photons', 'lidar_ratio']
        # Prepare for radiances
        if keyword == 'umu':
            self.n_umu = len(tokens)  # The number of umu values
            self.umu = np.array(tokens).astype(np.float)
        if keyword == 'phi':
            self.n_phi = len(tokens)
            self.phi = np.array(tokens).astype(np.float)
        if keyword == 'output_user':
            self.output_user = [token.replace('lambda', 'wvl') for token in tokens]  # lambda is a keyword
            self.output_user = [token.replace('wavenumber', 'wvn') for token in self.output_user]  # abbreviate wavenumber
            self.output_user = [token.replace('wavelength', 'wvl') for token in self.output_user]  # abbreviate wavelength
            self.fluxline = []
        if keyword == 'polradtran' and tokens[0] == 'nstokes':
            self.nstokes = int(tokens[1])
            self.prepare_for_polradtran()
        if keyword == 'zout' or keyword == 'zout_sea':  # Determine number of output levels
            if all([re.match(Case.re_isSingleOutputLevel, token.lower()) for token in tokens]):
                self.n_zout = len(tokens)
            else:  # Number of output levels is not known, will have to auto-detect
                self.n_zout = '?'

    def prepare_for_polradtran(self):
        """ Prepare for output from the polradtran solver
        :return:
        """
        stokescomps = ['I', 'Q', 'U', 'V']
        self.fluxline = ['wvl']
        for istokes in range(self.nstokes):
            self.fluxline.extend(['down_flux' + stokescomps[istokes], 'up_flux' + stokescomps[istokes]])

    def append_option(self, option, origin=('user', None)):
        ''' Append a libRadtran/uvspec options to this uvspec case. It will be appended at the end of the file
        :param option: A list containing the keyword and keyword parameters (tokens)
        :param origin: A 2-tuple giving the origin of the option and a "line number" reference. Default ('user', None)
        uvspec options.
        :return:
        '''
        self.options.append(option[0])  # the option keyword (string)
        self.tokens.append(option[1:])  # The tokens following the keyword (list of strings)
        self.filorigin.append(origin)  # The origin of this keyword
        self.optionobj.append(uvsOptions[option[0]])  # The option object
        # Make any possible preparations for occurance of this keyword
        self.prepare_for(option[0],option[1:])

    def alter_option(self, option, origin=('user', None)):
        """ Alter the parameters of a uvspec input option. If the option is not found, the option is appended with
        append_option instead.

        :param option: List of keyword and tokens (parameters) to provide to the option keyword (list of strings).
        :param origin: A 2-tuple noting the "origin" of the change to this keyword. Default ('user', None)
        :return:
        """
        ioption = self.options.index(option[0])
        if ioption:
            self.tokens[ioption] = option[1:]  # The tokens following the keyword (list of strings)
            self.filorigin[ioption] = origin  # The origin of this keyword
            self.prepare_for(option[0], option[1:])
        else:
            self.append_option(option, origin)  # just append the option if not found

    def del_option(self, option, all=True):
        """ Delete a uvspec input option matching the given option.

        :param option: Keyword of option to be deleted
        :param all: A flag indicating if all matching options must be deleted or only the first occurrence. The
        default is to delete all matching occurrences.
        :return: True if an option was deleted or False if not
        """

        #TODO consider providing warning if options does not exist
        deletedsomething = False
        while option in self.options:
            deletedsomething = True
            ioption = self.options.index(option)
            self.options.pop(ioption)
            self.tokens.pop(ioption)
            self.filorigin.pop(ioption)
            self.optionobj.pop(ioption)
        if deletedsomething:
            # Run through all options and reconstruct preparations
            for (ioption, option) in enumerate(self.options):
                self.prepare_for(self.options[ioption], self.tokens[ioption])
        return deletedsomething

    @staticmethod
    def read(path, includes_seen=[]):
        """ Reads a libRadtran input file. This will construct the libRadtran case from the contents of the .INP file
        Adapted from code by libRadtran developers.

        :param path: File path from which to read the uvspec input
        :param includes_seen: List of files already included (for recursion purposes to avoid infinite include loops)
        :return: data, line_nos, path
          where data is the full data in the file with includes, line_nos shows the source of every line and
          path is the path to the main input file.

        """
        # Get the full path
        path = os.path.abspath(path)
        folder = os.path.dirname(path)
        # print includes_seen
        with open(path, 'rt') as INPfile:
            opdata = INPfile.readlines()  # read in the entire file and process in memory afterwards
        opdata = [line.split() for line in opdata]  # Split lines into keywords and option parameters (tokens)
        line_nos = [(path, i) for i in xrange(1, len(opdata)+1)]
        # Remove lines with comments and merge continuous lines
        line = 0
        while line < len(opdata):  # Remove empty lines, comments and merge options split over multiple file lines
            # Skip empty lines
            if not opdata[line]:
                opdata.pop(line)
                line_nos.pop(line)
                continue
            # Skip comments  #TODO save comments as well in the librad.Case
            if opdata[line][0].startswith("#"):
                opdata.pop(line)
                line_nos.pop(line)
                continue
            # Remove comments from the line  #TODO save comments into the librad.Case
            elif [True for word in opdata[line] if (word.find("#") != -1)]:
                tmp = []
                for word in opdata[line]:
                    pos = word.find("#")
                    if pos != -1:
                        if word != "#":
                            tmp.append(word[:pos])
                        break
                    else:
                        tmp.append(word)
                opdata[line] = tmp
            # continuous line
            elif opdata[line][-1].endswith("\\"):
                opdata[line][-1] = opdata[pos][-1][:-1] # remove the \
                    # if the \ was preceded and continued by whitespace
                if opdata[line][-1] == "":
                    opdata[line].pop()
                    line_nos.pop()
                opdata[line].extend(opdata[line + 1])
            else:
                line += 1
        # Get the includes and include them at the point where the include keyword appears
        all_opdata = []
        all_line_nos = []
        this_includes = includes_seen # These are the include files seen up to this point in the recursion
        for line in xrange(len(opdata)):  # Run through all lines in the file and perform include substitutions
            if opdata[line][0] == 'include':  # Have found an include file, read it and append to all_opdata
                include_file = opdata[line][1]  # Obtain the filename
                include_file = os.path.normpath(os.path.join(folder, include_file))  # Extend to full path name
                if this_includes.count(include_file):  # This file has been included before
                    raise IOError('Attempted to include the same uvspec input file more than once. The included files are ', this_includes)
                this_includes.append(include_file)  # Add the filename to the list of included files already seen
                # Read the data from the included file
                inc_opdata, inc_line_nos, x = Case.read(include_file, this_includes)  # Recursion
                # Insert the data at this point in the file option and line number lists
                all_opdata.extend(inc_opdata)
                all_line_nos.extend(inc_line_nos)
            else:  # Not an include option, just append
                all_opdata.append(opdata[line])
                all_line_nos.append(line_nos[line])
        return all_opdata, all_line_nos, path

        # The following is the old code that inserted all include files at the end, which is not the same
        # behaviour as a C #INCLUDE statement, which inserts at the point where the #INCLUDE appears.
        # buff = 0
        # this_includes = []
        # for line in xrange(len(opdata)):
        #     if (opdata[line + buff][0] == "include" and len(opdata[line + buff]) == 2):
        #         include_file = opdata.pop(line + buff)[1]  # Obtain the filename of the included file and pop option
        #         include_file = os.path.normpath(os.path.join(folder, include_file))  # Extend to full path name
        #         this_includes.append(include_file)  # Add the filename to the list of included files at this level
        #         line_nos.pop(line + buff)  # Also pop the line numbers from the list of line numbers
        #         buff -= 1  # Take into account that a line has been removed
        # for include_path in this_includes:
        #     if not os.path.exists(include_path):
        #         msg = "Include file '%s' in '%s' does not exist." % (include_path, path)
        #         raise IOError(msg)
        #         #print " " * len(includes) + msg
        #     # If the file has been included before.
        #     elif (this_includes + includes_seen).count(include_path) != 1:  # Count number of times file is included
        #         msg = "File %s included more than once in %s. Please fix this." % (include_path, path)
        #         raise IOError(msg)
        #         #self.error_txt.append(msg)
        #         #print " " * len(includes) + msg
        #     else:
        #         include_data = Case.read(include_path, includes_seen + this_includes)
        #         # The include file might contain errors and return None.
        #         if include_data:
        #             opdata.extend(include_data[0])  #TODO surely an insert rather than extend ?
        #             line_nos.extend(include_data[1])
        # print opdata
        # print line_nos
        # return opdata, line_nos, path

    def __repr__(self):
        """ libRadtran/uvspec input data
        :return: The uvspec input data as it would appear in an input file.
        """
        uvsinp = []
        for (ioption, keyword) in enumerate(self.options):
            uvsinp.append(keyword + ' ' + ' '.join(self.tokens[ioption]))  #TODO add comments
        return '\n'.join(uvsinp)

    def write(self, filename=''):
        """ Write libRadtran/uvspec input to a file (.INP extension by default.
        If the filename input is given as '', a file save dialog will be presented
        :param filename: Filename to which to write the uvspec input
        :return:
        """
        if not filename:
            filename = easygui.filesavebox(msg='Please save the uvspec input file.', filetypes=["*.INP"])
        if filename[-4:].lower() != '.inp':
            filename += '.INP'
        with open(filename, 'wt') as uvINP:
            alldata = self.__repr__()
            uvINP.write(alldata)

    def distribute_flux_data(self, fluxdata):
        """ Distribute flux data read from uvspec output file to various fields

        :param fluxdata: Flux (irradiance) data read from uvspec output file
        :return:
        """
        # First split the data amongst output levels or output wavelengths/wavenumbers
        # We assume and handle only those cases of output_user where the primary variable is
        # zout, lambda (wvl) or wavenumber (wvn)

        if self.output_user:  # distribute to user-defined output variables
            fields = self.output_user
        else:
            fields = self.fluxline
        datashape = fluxdata.shape
        if len(datashape) == 1:  # Need some special cases to deal with single line output files
            linecount = 1
        else:
            linecount = datashape[1]
        if len(datashape) == 2:
            if datashape[1] != len(fields):  # number of fields in data does not match number of fields
                print datashape[1], ' not the same as ', len(fields) #TODO try to deal with this
        if fields[0] == 'zout':  # output level is the primary variable
            if self.n_zout == '?':  # Unknown number of output levels
                # Try just using the number of unique values in the first column
                self.n_zout = len(np.unique(fluxdata[:,0]))
            fluxdata = fluxdata.reshape((self.n_zout, -1, linecount), order='F')
        elif fields[0] == 'wvl' or fields[0] == 'wvn':  # wavelength/wavenumber is the primary variable
            if self.n_zout == '?':  # Don't know number of output levels
                self.n_wvl = len(np.unique(fluxdata[:,0]))  # Try to determine number of wavelengths/wavenumbers
                fluxdata = fluxdata.reshape((self.n_wvl, -1, linecount), order='F')
            else:
                fluxdata = fluxdata.reshape((-1, self.n_zout, linecount), order='F')
        else:  # Assume secondary variable is zout
            fluxdata = fluxdata.reshape((-1, self.n_zout, linecount), order='F')  #TODO provide warning or something
        self.fluxdata = fluxdata  # retain the flux data in the instance
        if linecount == 1:
            for (ifield, field) in enumerate(fields):
                setattr(self, field, np.squeeze(fluxdata[ifield]))
        else:
            for (ifield, field) in enumerate(fields):
                setattr(self, field, np.squeeze(fluxdata[:, :, ifield]))
        # Clean up zout and wavelength/wavenumber data
        self.zout = np.unique(self.zout)
        if len(self.zout) > 0:
            self.n_zout = len(self.zout)
        self.wvn = np.unique(self.wvn)
        self.wvl = np.unique(self.wvl)
        if len(self.wvl) > 0 and len(self.wvn) == 0:  # Calculate wavenumbers if wavelengths available
            self.wvn = 1e7 / self.wvl
        if len(self.wvn) > 0 and len(self.wvl) == 0:  # Calculate wavelengths if wavenumbers available
            self.wvl = 1e7 / self.wvn
        self.n_wvl = len(self.wvl)

    def readout(self, filename=None):
        """ Read uvspec output. The result is placed into a dictionary called self.out

         The general process of reading is:
          1) If the user has specified output_user, just assume a flat file and read using
             np.loadtxt or np.genfromtxt.
          2) If not user_output and the solver has no radiance blocks, assume a flat file and
             read using np.loadtxt. The variables to be read should be contained in the self.fluxline attribute.
          3) Otherwise, if the output has radiance blocks, read those depending on the radiance block
             format for the specific solver. Keep reading flux and radiance blocks until the file is exhausted.

         Once the data has all been read, the data is split up between the number of output levels and number
         of wavelengths. For radiance data, the order of numpy dimensions is umu, phi, wavelength and zout. That is,
         if a case has multiple zenith angles, multiple azimuth angles, multiple wavelengths and multiple output
         levels, the radiance property uu will have 4 dimensions.

        Output from uvspec depends on the solver and a number of other
        inputs, including the directive 'output_user'.
        For the solvers disort, sdisort, spsdisort and presumably also
        disort2, the irradiance (flux) outputs default to
          lambda edir edn eup uavgdir uavgdn uavgup

        If radiances (intensities) have been requested with the umu
        (cosine zenith angles input), each line of flux data is followed
        by a block of radiance data as follows:

         umu(0) u0u(umu(0))
         umu(1) u0u(umu(1))
         . . . .
         . . . .
         umu(n) u0u(umu(n))

        u0u is the azimuthally averaged radiance for the requested zenith
        angles.

        If azimuth angles (phi) have also been specified, then the
        radiance block is extended as follows:

                                phi(0)        ...     phi(m)
         umu(0) u0u(umu(0)) uu(umu(0),phi(0)) ... uu(umu(0),phi(m))
         umu(1) u0u(umu(1)) uu(umu(1),phi(0)) ... uu(umu(1),phi(m))
         . . . .
         . . . .
         umu(n) u0u(umu(n)) uu(umu(n),phi(0)) ... uu(umu(n),phi(m))

        Radiance outputs are not affected by output_user options.


        For the polradtran solver, the flux block is as follows:
           lambda down_flux(1) up_flux(1) ... down_flux(iS) up_flux(iS)

        where iS is the number of Stokes parameters specified using the
        'polradtran_nstokes' directive.
        If umu and phi are also specified, the radiance block is as
        follows:
                                 phi(0)      ...      phi(m)
        Stokes vector I
        umu(0) u0u(umu(0)) uu(umu(0),phi(0)) ... uu(umu(0),phi(m))
        umu(1) u0u(umu(1)) uu(umu(1),phi(0)) ... uu(umu(1),phi(m))
        . . . .
        . . . .
        umu(n) u0u(umu(n)) uu(umu(n),phi(0)) ... uu(umu(n),phi(m))
        Stokes vector Q
        . . .
        . . .
        u0u (azimuthally averaged radiance) is always zero for
        polradtran.

        For the two-stream solver (twostr), the flux block is
           lambda edir edn eup uavg

        The directive keyword 'brightness' can also change output. The
        documentation simply states that radiances and irradiances are
        just converted to brightness temperatures.

        The keyword directive 'zout' and it's parameters will influence
        output format as well. In general the output is repeated for each
        given value of zout or zout_sea.

        The keyword directive 'output' and its parameters will also have
        a major effect.

        'output sum'

        The keyword directive 'header' should not be used at all. This
        produces some header information in the output that will cause
        errors. An error is issued of the 'header' keyword is used in the
        input.

        :param filename: File from which to read the output. Defaults to name of input file, but with the .OUT
        extension.
        :return:
        """
        if filename is None:
            filename = self.outfile
        elif filename == '':
            filename = easygui.fileopenbox(msg='Please select the uvspec output file.', filetypes=["*.OUT"])
        fluxdata = []
        if self.output_user:
            fluxdata = np.loadtxt(filename)
            self.distribute_flux_data(fluxdata)
        elif (self.n_phi == 0 and self.n_umu == 0) or self.solver == 'sslidar':  # There are no radiance blocks
            fluxdata = np.loadtxt(filename)
            self.distribute_flux_data(fluxdata)
        # elif self.n_phi == 0:   # Not sure about format for n_umu > 0, n_phi == 0
        #  Look at example UVSPEC_FILTER_SOLAR.INP, which indicates manual is not correct
        #     fluxdata = np.loadtext(filename)
        #     # Take away the radiance data
        #     raddata = fluxdata[len(self.fluxline):]
        #     fluxdata = fluxdata[:len(self.fluxline)]
        #     self.raddata = raddata
        #     self.fluxdata = fluxdata
        #     self.distribute_flux_data(fluxdata)
        else:  # radiance blocks  #TODO polradtran has different format
            phicheck = []
            rad3D = []  # full 3D radiance data is here (umu, phi and wavelength)
            with open(filename, 'rt') as uvOUT:
                # Read and append a line of flux data
                txtline = uvOUT.readline();
                while txtline:
                    fluxline = np.array(txtline.split()).astype(np.float)
                    if len(fluxdata) > 0:
                        fluxdata = np.vstack((fluxdata, fluxline))
                    else:
                        fluxdata = fluxline
                    # Read the radiance block
                    if self.n_phi > 0:  # There is a line of phi angles
                        philine = uvOUT.readline()  #TODO check that phi angles are correct
                        phicheck = np.array(philine.split()).astype(np.float)
                    # Read the lines for the umu values
                    raddata = []
                    for i_umuline in range(self.n_umu):  #TODO this is wrong if there is no phi specified - see manual
                        umuline = uvOUT.readline()
                        radline = np.array(umuline.split()).astype(np.float)
                        if len(raddata) > 0:
                            raddata = np.vstack((raddata, radline))
                        else:
                            raddata = radline
                    if len(rad3D) > 0:
                        rad3D = np.dstack((rad3D, raddata))
                    else:
                        rad3D = raddata
                    txtline = uvOUT.readline()  # Read what should be the next line of flux data
            self.distribute_flux_data(fluxdata)  # distribute the flux data, which should also determine
                                                 # the number of wavelengths definitively
            self.rad3D = rad3D  # all the data
            self.phi_check = phicheck
            self.u0u = rad3D[:,1]
            #if self.u0u.shape[1] / self.n_wvl > 1:  #TODO problem here with single wavelength data
            #    self.u0u = self.u0u.reshape((self.n_umu, self.n_wvl, -1), order='F')
            if rad3D.ndim == 3:
                self.uu = rad3D[:,2:,:]
                # If there are multiple zout levels, must reshape self.uu
                if self.uu.shape[2] / self.n_wvl > 1:
                    self.uu = self.uu.reshape((self.n_umu, self.n_phi, self.n_wvl, -1), order='F')
            else:
                self.uu = rad3D[:,2:]
                if self.uu.shape[1] / self.n_wvl > 1:  # if multiple output levels, reshape the radiance data appropriately
                    self.uu = self.uu.reshape((self.n_umu, self.n_wvl, -1), order='F')



