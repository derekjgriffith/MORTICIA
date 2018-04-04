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

import warnings
import mitsuba.core as mitcor
import mitsuba.render as mitren

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
        v1 = mitcor.Vector4(x4x4[0, 0], x4x4[0, 1], x4x4[0, 2], x4x4[0, 3])
        v2 = mitcor.Vector4(x4x4[1, 0], x4x4[1, 1], x4x4[1, 2], x4x4[1, 3])
        v3 = mitcor.Vector4(x4x4[2, 0], x4x4[2, 1], x4x4[2, 2], x4x4[2, 3])
        v4 = mitcor.Vector4(x4x4[3, 0], x4x4[3, 1], x4x4[3, 2], x4x4[3, 3])
        xform = mitcor.Transform(mitcor.Matrix4x4(v1, v2, v3, v4))
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
        origin_position = mitcor.Point(origin_position[0],
                                             origin_position[1],
                                             origin_position[2])
        target_position = mitcor.Point(target_position[0],
                                             target_position[1],
                                             target_position[2])
        up_direction = mitcor.Vector(up_direction[0],
                                           up_direction[1],
                                           up_direction[2])
        self.xform = mitcor.Transform.lookAt(origin_position, target_position, up_direction) * self.xform

    def translate(self, translation):
        """

        :param translation:
        :return:
        """
        translation = mitcor.Vector3(translation[0],
                                           translation[1],
                                           translation[2])
        translation = mitcor.Transform.translate(translation)
        self.xform = translation * self.xform

    def rotate(self, axis=(0.0, 0.0, 1.0), angle=0.0):
        axis = mitcor.Vector3(axis[0], axis[1], axis[2])
        self.xform = mitcor.Transform.rotate(axis, angle) * self.xform

    def scale(self, scale_factors=(1.0, 1.0, 1.0)):
        if len(scale_factors) == 3:
            scale_factors = mitcor.Vector3(scale_factors[0],
                                                 scale_factors[1],
                                                 scale_factors[2])
        else:
            scale_factors = mitcor.Vector3(scale_factors,
                                                 scale_factors,
                                                 scale_factors)

        self.xform = mitcor.Transform.scale(scale_factors) * self.xform

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
        self.fileResolver = mitcor.Thread.getThread().getFileResolver()
        if scenefolder is not None:  # the user input something
            if scenefolder:  #  not blank
                if isinstance(scenefolder, basestring):
                    self.fileResolver.appendPath(scenefolder)
                else:  # assume a list of strings
                    for path in scenefolder:
                        self.fileResolver.appendPath(path)
        theparamMap = mitcor.StringMap()
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
            self.scene = mitren.SceneHandler.loadScene(
                           self.fileResolver.resolve(scenefile), theparamMap)
        else: # Create a new scene
            self.scene = mitren.Scene()
        # Add a plugin manager
        self.plugin_mngr = mitcor.PluginManager.getInstance()

    def add_directional_light(self, direction=(0.0, 0.0, -1.0), irradiance=1.0, samplingWeight=1.0):
        """
        Add a directional light (emitter) to a Mitsuba scene instance. Every call will add another emitter.
        :param direction: This is a 3-vector providing the direction of the light. This can
        be a list or tuple of floats, or a numpy vector. Defaults to (0.0, 0.0, -1.0), that is with
        the light pointing in the negative z direction (towards nadir, as for the sun in the zenith).
        :param irradiance: A scalar or vector of numbers providing the irradiance of the light source in canonical
        irradiance units. Defaults to 1.0.
        :param samplingWeight: Relative amount of samples allocated to this light source. Default is
        1.0.
        :return:
        """
        direction = mitcor.Vector3(direction[0], direction[1], direction[2])
        if isinstance(irradiance, np.ndarray):
            irradiance = irradiance.tolist()
        irradiance = mitcor.Spectrum(irradiance)
        dir_light_props = mitcor.Properties('directional')  # Carries properties of a directional emitter
        dir_light_props['direction'] = direction
        dir_light_props['irradiance'] = irradiance
        dir_light_props['samplingWeight'] = samplingWeight
        dir_light = self.plugin_mngr.createObject(dir_light_props)
        dir_light.configure()
        self.scene.addChild('directional', dir_light)

    def add_direct_sun(self, sza=0.0, saz=0.0, irradiance=1.0, samplingWeight=1.0):
        """
        Add the sun as a directional Mitsuba emitter, using solar zenith angle and azimuth angle.
        Every call will add another directional emitter to the Mitsuba scene
        :param sza: Solar zenith angle in degrees, defaults to 0.0 (sun in the zenith)
        :param saz: Solar azimuth angle in degrees, defaults to 0.0 (sun in the north). Measured from north through
        East. Mitsuba environments are defined with +Y to the zenith, but in `MORTICIA`, +Z is towards the zenith +X
        is North and +Y is East. Also note that the phi0 input for libRadtran is the azimuthal direction in which
        sunlight is moving.
        :param irradiance: A scalar or vector of numbers providing the irradiance of the light source in canonical
        irradiance units. Defaults to 1.0.
        :param samplingWeight: Relative amount of samples allocated to this light source. Default is
        1.0 and this is the value generally used in MORTICIA.
        :return:
        """
        # Just calculate the direction in the required coordinate frame add the directional light
        altitude = np.pi/2.0 - np.deg2rad(sza)
        azimuth = np.deg2rad(saz)
        z = np.sin(altitude)
        hyp = np.cos(altitude)
        y = hyp * np.cos(azimuth)
        x = hyp * np.sin(azimuth)
        light_direction = (-x, -y , -z)  # Light travels in opposite direction to vector towards the sun
        self.add_directional_light(direction=light_direction, irradiance=irradiance, samplingWeight=samplingWeight)

    def add_radiant_environment_map(self, filename, scale=1.0, toWorld=None, gamma=1.0, cache=None, samplingWeight=1.0):
        """
        Add a radiant enviroment map (envmap) to a Mitsuba scene taken from a file (typically OpenEXR format).
        :param filename: Filename of the file to fetch the REM from. Supports any image format supported by Mitsuba.
        For `MORTICIA` purposes, only use OpenEXR.
        :param scale: Scale the REM by this scalar value. Defaults to 1.0.
        :param toWorld: Transformation (only rotation matters) of the environment map. By default, the REM will be
        rotated so that the zenith is along the X-axis, the X-axis is towards the North and the Y-axis is towards the
        East.
        :param gamma: Override the gamma value of the REM source bitmap. For MORTICIA purposes, the REM should be
        absolute and high dynamic range (OpenEXR) and gamma will default to 1.0.
        :param cache: Set True to force MIP mapping of the REM and False to inhibit. Default is automatic which will
        cache the MIP map for images larger than 1 megapixel.
        :param samplingWeight: Relative weight to assign to this emitter. Default is 1.0. In MORTICIA we generally
        use the default.
        :return:
        """
        if self.scene.hasEnvironmentEmitter():
            warnings.warn('Mitsuba scene already has a radiant environment emitter. Only one is permitted.')
        if toWorld is None:  # Generate the default REM rotation, which is +Z towards the zenith
            toWorld = Transform()  # Get identity transform
            toWorld.rotate((1.0, 0.0, 0.0), 90.0)  # Rotate by 90 degrees about the x-axis
        envmap_props = mitcor.Properties('envmap')
        envmap_props['filename'] = filename
        envmap_props['toWorld'] = toWorld.xform
        envmap_props['gamma'] = gamma
        if cache is not None:
            envmap_props['cache'] = cache
        envmap_props['samplingWeight'] = samplingWeight
        envmap = self.plugin_mngr.createObject(envmap_props)
        envmap.configure()
        self.scene.addChild('environment', envmap)


