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

import writeLex  # This imports all the libradtran option definitions
import os
import easygui  # For file open dialogs
import numpy as np

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
        'lambda': 'Output wavelengths [nm]',
        'wavelength': 'Presumably also wavelengths [nm]',
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

    def __init__(self, casename='', filename=None, optiondict=None):
        """ A libRadtran/uvspec case
        :param casename: A user-defined name for the libRadtran/uvspec case
        :param filename: An optional filename from which to read the libRadtran/uvspec input
        :param optiondict: A dictionary of uvspec option names, with a list of text tokens for each option.
        e.g. {'sza': ['25.0']} will define the solar zenith angle option at 25 degrees. Options that are not found in
        the current dictionary will return errors. Besides that, no error checking is performed automatically.
        :return:
        """
        self.name = casename
        self.error_txt = []
        self.options = []  # options is a list [option_name (string), option_tokens (list of strings),
        self.tokens = []
        self.optionobj = []
        self.filorigin = []
        self.solver = 'disort'  # default, modified b y the rte_solver option keyword
        self.fluxline = ['lambda', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
        self.output_user = ''  # set with the output_user keyword
        self.output_quantity = '' # default is radiances and irradiances
        self.n_umu = 0
        self.n_phi = 0
        self.nstokes = 1  # default number of stokes parameters for polradtran
        if filename is not None:
            if filename[-4:].lower() != '.inp':
                filename = filename + '.INP'
            data, line_nos, path = Case.read(path=filename)
            self.infile = path
            self.outfile = self.infile[:-4] + '.OUT'
            # option_object, source_file_nos (filename, line_number)]
            # Process the results into lists of options
            for (optnumber, option) in enumerate(data):
                self.options.append(option[0])  # the option keyword (string)
                self.tokens.append(option[1:])  # The tokens following the keyword (list of strings)
                self.filorigin.append(line_nos[optnumber])  # The file origin of this keyword
                self.optionobj.append(uvsOptions[self.options[optnumber]])  # The option object
                # Make any possible preparations for occurance of this keyword
                self.prepare_for(option[0],option[1:])


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
            self.fluxline = ['lambda', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
        elif any([self.solver == thesolver for thesolver in ['twostr','rodents']]):
            self.fluxline = ['lambda', 'edir', 'edn', 'eup', 'uavg']
        elif self.solver == 'polradtran':
            self.prepare_for_polradtran()
        elif self.solver == 'sslidar':
            self.fluxline = ['center_of_range', 'number_of_photons', 'lidar-ratio']
        # Prepare for radiances
        if keyword == 'umu':
            self.n_umu = len(tokens)  # The number of umu values
        if keyword == 'phi':
            self.n_phi = len(tokens)
        if keyword == 'output_user':
            self.output_user = tokens
            self.fluxline = []
        if keyword == 'polradtran' and tokens[0] == 'nstokes':
            self.nstokes = int(tokens[1])
            self.prepare_for_polradtran()

    def prepare_for_polradtran(self):
        """ Prepare for output from the polradtran solver
        :return:
        """
        stokescomps = ['I', 'Q', 'U', 'V']
        self.fluxline = ['lambda']
        for istokes in range(self.nstokes):
            self.fluxline.extend(['down_flux' + stokescomps[istokes], 'up_flux' + stokescomps[istokes]])

    def addoption(self, optiondict):
        ''' Add libRadtran/uvspec options to this uvspec case.
        :param optiondict: A dictionary of options to add to the case. The dict keys must be valid names of libRadtran/
        uvspec options.
        :return:
        '''

    @staticmethod
    def read(path='', includes=[]):
        """ Reads a libRadtran input file. This will construct the libRadtran case from the contents of the .INP file
        Adapted from code by libRadtran developers.

        """
        if not path:
            # Open a dialog
            path = easygui.fileopenbox(msg='Please select a uvspec input file.', filetypes=["*.INP"])
        with open(path, 'rt') as INPfile:
            data = INPfile.readlines()
        data = [line.split() for line in data]
        line_nos = [(path, i) for i in xrange(1, len(data)+1)]
        # Remove lines with comments and merge continuous lines
        line = 0
        while line < len(data):
            # Skip empty lines
            if not data[line]:
                data.pop(line)
                line_nos.pop(line)
                continue
            # Skip comments
            if data[line][0].startswith("#"):
                data.pop(line)
                line_nos.pop(line)
                continue
            # Remove comments from the line  #TODO save comments into the librad.Case
            elif [True for word in data[line] if (word.find("#") != -1)]:
                tmp = []
                for word in data[line]:
                    pos = word.find("#")
                    if pos != -1:
                        if word != "#":
                            tmp.append(word[:pos])
                        break
                    else:
                        tmp.append(word)
                data[line] = tmp
            # continous line
            elif data[line][-1].endswith("\\"):
                data[line][-1] = data[pos][-1][:-1] # remove the \
                    # if the \ was preceded and continued by whitespace
                if data[line][-1] == "":
                    data[line].pop()
                    line_nos.pop()
                data[line].extend(data[line + 1])
            else:
                line += 1
        # Get the includes and include them
        buff = 0
        this_includes = []
        for line in xrange(len(data)):
            if (data[line + buff][0] == "include" and len(data[line + buff]) == 2):
                this_includes.append(data.pop(line + buff)[1])
                line_nos.pop(line + buff)
                buff -= 1
        for include_path in this_includes:
            if not os.path.exists(path):
                msg = "Include file '%s' in '%s' does not exist." % (include_path, path)
                # self.error_txt.append(msg)
                #print " " * len(includes) + msg
            # If the file has been included before.
            elif (this_includes + includes).count(include_path) != 1:
                msg = "File %s included more than once in %s. Please fix this." % (include_path, path)
                #self.error_txt.append(msg)
                #print " " * len(includes) + msg
            else:
                include_data = Case.read(include_path, includes + this_includes)
                # The include file might contain errors and return None.
                if include_data:
                    data.extend(include_data[0])
                    line_nos.extend(include_data[1])
        return data, line_nos, path

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

    def determine_out_fmt(self):
        """ Runs through all input options and determines the output format
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

        The keyword diretive 'header' should not be used at all. This
        produces some header information in the output that will cause
        errors. An error is issued of the 'header' keyword is used in the
        input.
        :return:
        """
        pass



    def readout(self, filename=None):
        """ Read uvspec output. The result is placed into a dictionary called self.out

        :param filename: File from which to read the output. Defaults to name of input file, but with the .OUT
        extension.
        :return:
        """
        if filename is None:
            filename = self.outfile
        elif filename == '':
            filename = easygui.fileopenbox(msg='Please select the uvspec output file.', filetypes=["*.OUT"])
        if self.output_user:
            data = np.loadtxt(filename)  #TODO sutff here
            self.data = data
        elif:  #TODO
