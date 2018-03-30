__author__ = 'DGriffith'
"""
 *--------------------------------------------------------------------
 *
 * This file is part of MORTICIA.
 * Copyright (c) 2015-2018 by Derek Griffith
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
Class for Mitsuba support
"""
#import lxml
import numpy as np
# The following is messy, but nothing else seems to work.
# Putting Mitsuba path into morticia.pth does not work.
import os, sys
# NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
# '\' to avoid pitfalls with string escaping)
# Configure the search path for the Python extension module
sys.path.append('D:/Projects/MORTICIA/Render/Mitsuba 0.5.0/python/2.7')
# Ensure that Python will be able to find the Mitsuba core libraries
os.environ['PATH'] = 'D:/Projects/MORTICIA/Render/Mitsuba 0.5.0' + os.pathsep + os.environ['PATH']
import mitsuba

class Transform(object):
    """
    Encapsulates transformation objects for moving and rotating Mitsuba scene components
    """
    def __init__(self, x4x4=np.array([[1.0, 0.0, 0.0, 0.0],
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, 1.0, 0.0],
                                      [0.0, 0.0, 0.0, 1.0]])):
        """
        Construct a Mitsuba transformation object.
        :param x4x4: Numpy 4x4 array
        Defaults to the identity transformation.
        :return: Class Transform object with property xform, which contains the Mitsuba 4x4 matrix transform.
        """
        v1 = mitsuba.core.Vector4(x4x4[0, 0], x4x4[0, 1], x4x4[0, 2], x4x4[0, 3])
        v2 = mitsuba.core.Vector4(x4x4[1, 0], x4x4[1, 1], x4x4[1, 2], x4x4[1, 3])
        v3 = mitsuba.core.Vector4(x4x4[2, 0], x4x4[2, 1], x4x4[2, 2], x4x4[2, 3])
        v4 = mitsuba.core.Vector4(x4x4[3, 0], x4x4[3, 1], x4x4[3, 2], x4x4[3, 3])
        xform = mitsuba.core.Transform(mitsuba.core.Matrix4x4(v1, v2, v3, v4))
        self.xform = xform

    def toWorld_lookAt(self, origin_position=(10.0, 10.0, 10.0),
                         target_position=(0.0, 0.0, 0.0),
                         up_direction=(0.0, 0.0, 1.0)):
        """
        Construct a Mitsuba transformation object using origin and target positions in space together with an "up"
        direction to define e.g. camera position with up direction in the field of view. The current xform property is
        left multiplied by the new lookAt matrix.
        :param origin_position:
        :param target_position:
        :param up_direction:
        :return:
        """
        origin_position = mitsuba.core.Point(origin_position[0],
                                             origin_position[1],
                                             origin_position[2])
        target_position = mitsuba.core.Point(target_position[0],
                                             target_position[1],
                                             target_position[2])
        up_direction = mitsuba.core.Vector(up_direction[0],
                                           up_direction[1],
                                           up_direction[2])
        self.xform = mitsuba.core.Transform.lookAt(origin_position, target_position, up_direction) * self.xform

    def translate(self, translation):
        """

        :param translation:
        :return:
        """
        translation = mitsuba.core.Vector3(translation[0],
                                           translation[1],
                                           translation[2])
        translation = mitsuba.core.Transform.translate(translation)
        self.xform = translation * self.xform

    def rotate(self, axis=(0.0, 0.0, 1.0), angle=0.0):
        axis = mitsuba.core.Vector3(axis[0], axis[1], axis[2])
        self.xform = mitsuba.core.Transform.rotate(axis, angle) * self.xform

    def scale(self, scale_factors=(1.0, 1.0, 1.0)):
        if len(scale_factors) == 3:
            scale_factors = mitsuba.core.Vector3(scale_factors[0],
                                                 scale_factors[1],
                                                 scale_factors[2])
        else:
            scale_factors = mitsuba.core.Vector3(scale_factors,
                                                 scale_factors,
                                                 scale_factors)

        self.xform = mitsuba.core.Transform.scale(scale_factors) * self.xform

    def __mul__(self, other):
        self.xform = self * other

    def __rmul__(self, other):
        self.xform = other * self


class Mitsuba(object):
    """ The Mitsuba class encapsulates the information that can be contained in a Mitsuba scene file.
    This class contains the information in the scene file in an accessible format for the purpose
    of manipulating Mitsuba scenes, writing the updated scene file and performing rendering
    in a multiprocessing environment.

    """
    def __init__(self, scenefile=None, paramMap=None, scenefolder=None):
        """

        :param scenefile: Name of file from which to load a scene (string). Must not include the path, but must have
        the extension .xml.
        :param paramMap: Dictionary of parameters which will get substituted when the scene is loaded.
        :param scenefolder: Path to the folder containing the scene if not in current directory. A list of folders in
        which to search for the scene can also be provided (string or list of strings).
        :return:
        """

        # Set up the file resolver for this Mitsuba object
        # Get a reference to the threads file resolver
        self.fileResolver = mitsuba.core.Thread.getThread().getFileResolver()
        if scenefolder is not None:  # the user input something
            if scenefolder:  #  not blank
                if isinstance(scenefolder, basestring):
                    self.fileResolver.appendPath(scenefolder)
                else:  # assume a list of strings
                    for path in scenefolder:
                        self.fileResolver.appendPath(path)
        theparamMap = mitsuba.core.StringMap()
        if paramMap is not None:  # dict provided, convert to StringMap
            for param, value in paramMap.iteritems():
                theparamMap[param] = value
        # Read in the scene file if requested
        if scenefile is not None:
            if not scenefile:
                # Open a dialog to get the filename
                import easygui
                scenefile = easygui.fileopenbox(msg='Please select a Mitsuba scene file.', filetypes=["*.xml"])
            # Load the scene
            self.scene = mitsuba.render.SceneHandler.loadScene(
                           self.fileResolver.resolve(scenefile), theparamMap)

    def add_directional_light(self, tranform, direction, irradiance, samplingWeight=1.0):
        pass

