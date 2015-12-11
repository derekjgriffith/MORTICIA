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
There is no Windows-native version of libRadtran, so that generally implies that you are running on a Unix/Linux
machine or a virtual Unix-like environment on Windows (e.g. Cygwin).
The code here is taken from libRadtran version 2.0

"""

import writeLex  # This imports all the libradtran option definitions
import os
import easygui  # For file open dialogs

class Case():
    ''' Class which encapsulates a run case of libRadtran/uvspec.
    This class has methods to read libRadtran/uvspec input files, write uvspec input files, run uvspec in parallel on
    multiple compute nodes and read uvspec output files. An important use-case is that of reading a uvspec input
    file called the "base case", altering the parameters of particular option keywords and then running the case
    and reading the outputs. This class is also used by the RadEnv class which encapsulates a radiant environment.
    Construction of radiant environment maps typically requires running an array of librad.Case instances.
    '''
    def __init__(self, filename=None, optiondict={}):
        ''' A libRadtran/uvspec case
        :param filename: An optional filename from which to read the libRadtran/uvspec input
        :param optiondict: A dictionary of uvspec option names, with a list of text tokens for each option.
        e.g. {'sza': ['25.0']} will define the solar zenith angle option at 25 degrees. Options that are not found in
        the current dictionary will return errors. Besides that, no error checking is performed automatically.
        :return:
        '''
        self.error_txt = []
        if filename is not None:
            data, line_nos = Case.read(path=filename)
            self.data = data
            self.line_nos = line_nos

    def addoption(self, optiondict):
        ''' Add libRadtran/uvspec options to this uvspec case.
        :param optiondict: A dictionary of options to add to the case. The dict keys must be valid names of libRadtran/
        uvspec options.
        :return:
        '''

    @staticmethod
    def read(path='', includes=[]):
        ''' Reads a libRadtran input file. This will construct the libRadtran case from the contents of the .INP file
        Adapted from code by libRadtran developers.

        '''
        if not path:
            # Open a dialog
            path = easygui.fileopenbox(msg='Please select a uvspec input file.', filetypes=["*.INP"])
        # try:
        #     f = open(path, "r")
        #     data = f.readlines()
        #     f.close()
        # except IOError:
        #     msg = "Could not open '%s'." % (path)
        #     #self.error_txt.append(msg)
        #
        #     return
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
            # Remove comments from the line
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
            if (data[line + buff][0] == "include" and \
                    len(data[line + buff]) == 2):
                this_includes.append(data.pop(line + buff)[1])
                line_nos.pop(line + buff)
                buff -= 1
        for include_path in this_includes:
            if not os.path.exists(path):
                msg = "Include '%s' in '%s' does not exist." % (include_path, path)
                # self.error_txt.append(msg)
                #print " " * len(includes) + msg
            # If the file has been included before.
            elif (this_includes + includes).count(include_path) != 1:
                msg = "Include %s included more than one time." + \
                    " One time in %s. Please fix this " % (include_path, path)
                #self.error_txt.append(msg)
                #print " " * len(includes) + msg
            else:
                include_data = Case.read(include_path,
                                             includes + this_includes)
                # The include file might contains errors and return None.
                if include_data:
                    data.extend(include_data[0])
                    line_nos.extend(include_data[1])

        # Only Display error message in the most outer scope.
        #if self.error_txt and not includes:
        #    ErrorMessage(self.error_txt)
        #    self.error_txt = []
        #self.data = data
        #self.line_nos = line_nos
        return data, line_nos

