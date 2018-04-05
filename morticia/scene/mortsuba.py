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
import warnings
import mitsuba.core as mitcor
import mitsuba.render as mitren
import re
# Read the Mitsuba manual section on Python integration. On windows it is necessary to explicitly specify the location
# of the Mitsuba installation as follows:

import os, sys
if sys.platform == 'win32':
    # NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
    # '\' to avoid pitfalls with string escaping)
    # Configure the search path for the Python extension module
    # Replace strings below with location of your Mitsuba
    sys.path.append('D:/Projects/MORTICIA/Render/Mitsuba 0.5.0/python/2.7')
    # Ensure that Python will be able to find the Mitsuba core libraries
    os.environ['PATH'] = 'D:/Projects/MORTICIA/Render/Mitsuba 0.5.0' + os.pathsep + os.environ['PATH']
else:
    pass  # Assume that user has done the necessary
    # On Linux, ensure that the `source setpath.sh` appears in your .bashrc


def probe_mitsuba_SPECTRUM_SAMPLES():
    """
    Determine the number of SPECTRUM_SAMPLES for which this Mitsuba has been compiled.
    This is a messy and indirect solution and the Python API may have a direct method hidden somewhere.
    :return: The number of SPECTRUM_SAMPLES used when the available Mitsuba binary was compiled. See the Mitsuba
    manual for further details on the SPECTRUM_SAMPLES compilation switch.
    """
    # Instantiate a spectrum with 2 samples and process the error message
    spectrum_samples = None
    scan_match = re.compile('main \[core\.cpp:\d*\] Spectrum: expected (\d*) arguments')
    try:
        myspectrum = mitcor.Spectrum([1.0, 1.0])  # 2-sample spectrum should be rejected with Mitsuba error message
    except RuntimeError, error:
        spectrum_samples = int(re.match(scan_match, error.message).groups()[0])
    return spectrum_samples

# Set up some defaults
SPECTRUM_SAMPLES = probe_mitsuba_SPECTRUM_SAMPLES()
if SPECTRUM_SAMPLES == 3:
    default_pixelFormat = 'rgb'
else:
    default_pixelFormat = 'spectrum'

default_banner = False  # Do not, by default put Mitsuba banner on output images
default_componentFormat = 'float32'  # Use maximum accuracy for output Mitsuba image formats which allow this
default_attachLog = False  # Do not, by default attach complete rendering log to output hdrfilm images
default_height = 768  # Default film output height in pixels for Mitsuba renders
default_width = 576  # Default film output width in pixels for Mitusba renders
default_hdr_fileFormat = 'openexr'  # Default high dynamic range film output
default_highQualityEdges = True  # Use high quality edges by default in case image must be inserted into another
default_ldr_fileFormat = 'png'  # Default low dynamic range file format output
default_ldr_tonemapMethod = 'gamma'
default_ldr_gamma = -1.0
default_ldr_exposure = 0.0
default_ldr_key = 0.18
default_ldr_burn = 0.0

# Get a global Mitsuba plugin manager
plugin_mngr = mitcor.PluginManager.getInstance()

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

class Animation(object):
    """
    This class encapsulates animation transformations for Mitsuba objects (mainly target geometry components and
    sensor locations). Animations provide a series of transforms at a sequence of times. The time is assumed to be in
    seconds relative to the epoch (start time) of the scenario. The animation transformation is then a list of
    Transform objects with associated elapsed time in seconds.

    Unclear right now if Mitsuba Python bindings include the render Animation transformations.
    """
    pass

class ReconstructionFilter(object):

    def __init__(self, type='gaussian', stddev=2.0, B=1.0/3.0, C=1.0/3.0, lobes=3):
        self.type = type
        self.plugin_mngr = mitcor.PluginManager.getInstance()
        rfilter_props = mitcor.Properties(type)
           # Don't know what this property is or how to find out
        if type == 'mitchell' or type == 'catmullrom':
            rfilter_props['B'] = B
            rfilter_props['C'] = C
        if type == 'lanczos':
            rfilter_props['lobes'] = lobes
        self.rfilter = self.plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()

# Classes for different Mitsuba reconstruction filters
"""
The Filter classes encapsulate Mitsuba reconstruction filters applied to images after rendering. From the Mitsuba
manual:
Image reconstruction filters are responsible for converting a series of radiance samples generated
jointly by the sampler and integrator into the final output image that will be written to disk at the
end of a rendering process. This section gives a brief overview of the reconstruction filters that are
available in Mitsuba.There is no universally superior filter, and the final choice depends on a trade-off
between sharpness, ringing, and aliasing, and computational efficiency. Mitsuba currently has 6 reconstruction
filters, namely `box`, `tent`, `gaussian`, `mitchell`, `catmullrom` and `lanczos`.
"""
class BoxFilter(object):
    """
    The fastest, but also about the worst possible reconstruction filter, since it is extremely
    prone to aliasing. It is included mainly for completeness, though some rare situations
    may warrant its use. The filter radius defaults to 0.5 pixels.
    Note that the radius input to the BoxFilter constructor may not work.
    """
    def __init__(self, radius=0.5):
        self.filtertype = 'box'
        rfilter_props = mitcor.Properties(self.filtertype)
        rfilter_props['radius'] = radius
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()


class TentFilter(object):
    """
    Simple tent, or triangle filter.This reconstruction filter never suffers from ringing
    and usually causes less aliasing than a naive box filter.When rendering scenes with sharp brightness
    discontinuities, this may be useful; otherwise, negative-lobed filters will be preferable (e.g.
    Mitchell-Netravali or Lanczos Sinc). The filter radius defaults to 1.0 pixels.
    """
    def __init__(self):
        self.filtertype = 'tent'
        rfilter_props = mitcor.Properties(self.filtertype)
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()


class GaussianFilter(object):
    """
    This is a windowed Gaussian filter with configurable standard deviation.
    It produces pleasing results and never suffers from ringing, but may occasionally introduce too
    much blurring. When no reconstruction filter is explicitly requested, this is the default choice
    in Mitsuba. The gaussian standard deviation defaults to 0.5 pixels. Providing another value may not work.
    """
    def __init__(self, stddev=0.5):
        self.filtertype = 'gaussian'
        rfilter_props = mitcor.Properties(self.filtertype)
        rfilter_props['stddev'] = stddev
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()


class MitchellFilter(object):
    """
    Separable cubic spline reconstruction filter by Mitchell and
    Netravali [32]This is often a good compromise between sharpness and ringing.
    The plugin has two float-valued parameters named B and C that correspond to the two parameters
    in the original research paper. By default, these are set to the recommended value of
    1/3, but can be tweaked if desired.
    """
    def __init__(self, B=1.0/3.0, C=1.0/3.0):
        self.filtertype = 'mitchell'
        rfilter_props = mitcor.Properties(self.filtertype)
        rfilter_props['B'] = B
        rfilter_props['C'] = C
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()


class CatmullRomFilter(object):
    """
    This is a special version of the Mitchell-Netravali filter that has
    the constants B and C adjusted to produce higher sharpness at the cost of increased susceptibility
    to ringing.
    """
    def __init__(self):
        self.filtertype = 'catmullrom'
        rfilter_props = mitcor.Properties(self.filtertype)
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()


class LanczosFilter(object):
    """
    This is a windowed version of the theoretically optimal low-pass filter.
    It is generally one of the best available filters in terms of producing sharp high-quality
    output. Its main disadvantage is that it produces strong ringing around discontinuities, which
    can become a serious problem when rendering bright objects with sharp edges (for instance, a
    directly visible light source will have black fringing artifacts around it).This is also the computationally
    slowest reconstruction filter.
    This plugin has an integer-valued parameter named lobes, that sets the desired number of
    filter side-lobes.The higher, the closer the filter will approximate an optimal low-pass filter, but
    this also increases the susceptibility to ringing. Values of 2 or 3 are common (3 is the default).
    """
    def __init__(self, lobes=3):
        self.filtertype = 'lanczos'
        rfilter_props = mitcor.Properties(self.filtertype)
        rfilter_props['radius'] = radius
        rfilter_props['lobes'] = lobes
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()

# End of Mitsuba reconstruction filter classes

# Classes for various times of Mitsuba sensor "film".
"""
The Film classes encapsulates Mitsuba films. From the Mitsuba manual:
A film defines how conducted measurements are stored and converted into the final output file that
is written to disk at the end of the rendering process. Mitsuba comes with a few films that can write
to high and low dynamic range image formats (OpenEXR, JPEG or PNG), as well more scientifically
oriented data formats (e.g. Numpy, MATLAB or Mathematica).
"""

class HdrFilm(object):
    def __init__(self, width=default_width, height=default_height, cropOffsetX=None, cropOffsetY=None,
                 cropWidth=None, cropHeight=None, fileFormat=default_hdr_fileFormat, pixelFormat=default_pixelFormat,
                 componentFormat=default_componentFormat, attachLog = default_attachLog,
                 banner=default_banner, highQualityEdges=default_highQualityEdges, rfilter=None):
        """

        :param width:
        :param height:
        :param cropOffsetX:
        :param cropOffsetY:
        :param cropWidth:
        :param cropHeight:
        :param fileFormat:
        :param pixelFormat:
        :param componentFormat:
        :param attachLog:
        :param banner:
        :param highQualityEdges:
        :param rfilter:
        :return:
        """
        hdrfilm_props = mitcor.Properties('hdrfilm')
        hdrfilm_props['width'] = width
        hdrfilm_props['height'] = height
        hdrfilm_props['fileFormat'] = fileFormat
        if cropOffsetX is not None:
            hdrfilm_props['cropOffsetX'] = cropOffsetX
            hdrfilm_props['cropOffsetY'] = cropOffsetY
            hdrfilm_props['cropWidth'] = cropWidth
            hdrfilm_props['cropHeight'] = cropHeight
        hdrfilm_props['pixelFormat'] = pixelFormat
        hdrfilm_props['componentFormat'] = componentFormat
        hdrfilm_props['attachLog'] = attachLog
        hdrfilm_props['banner'] = banner
        hdrfilm_props['highQualityEdges'] = highQualityEdges
        hdrfilm = plugin_mngr.createObject(hdrfilm_props)
        if rfilter is not None:
            hdrfilm.addChild(rfilter)
        hdrfilm.configure()
        self.film = hdrfilm


class TiledHdrFilm(object):
    def __init__(self, width=default_width, height=default_height, cropOffsetX=None, cropOffsetY=None,
                 cropWidth=None, cropHeight=None, pixelFormat=default_pixelFormat,
                 componentFormat=default_componentFormat, rfilter=None):
        """

        :param width:
        :param height:
        :param cropOffsetX:
        :param cropOffsetY:
        :param cropWidth:
        :param cropHeight:
        :param pixelFormat:
        :param componentFormat:
        :param rfilter:
        :return:
        """
        hdrfilm_props = mitcor.Properties('tiledhdrfilm')
        hdrfilm_props['width'] = width
        hdrfilm_props['height'] = height
        if cropOffsetX is not None:
            hdrfilm_props['cropOffsetX'] = cropOffsetX
            hdrfilm_props['cropOffsetY'] = cropOffsetY
            hdrfilm_props['cropWidth'] = cropWidth
            hdrfilm_props['cropHeight'] = cropHeight
        hdrfilm_props['pixelFormat'] = pixelFormat
        hdrfilm_props['componentFormat'] = componentFormat
        hdrfilm = plugin_mngr.createObject(hdrfilm_props)
        if rfilter is not None:
            hdrfilm.addChild('filter', rfilter)
        hdrfilm.configure()
        self.film = hdrfilm

class LdrFilm(object):
    def __init__(self, width=default_width, height=default_height, cropOffsetX=None, cropOffsetY=None,
                 cropWidth=None, cropHeight=None, fileFormat=default_ldr_fileFormat,
                 pixelFormat=default_pixelFormat, tonemapMethod=default_ldr_tonemapMethod,
                 gamma=default_ldr_gamma, exposure=default_ldr_exposure, key=default_ldr_key,
                 burn=default_ldr_burn,
                 banner=default_banner, highQualityEdges=default_highQualityEdges, rfilter=None):
        """

        :param width: Width of sensor film in pixels
        :param height: Height of sensor film in pixels
        :param cropOffsetX: The crop parameters can optionally be provided to select a subrectangle
        of the output. In this case, Mitsuba will only render
        :param cropOffsetY:
        :param cropWidth:
        :param cropHeight:
        :param fileFormat: The desired output file format: png or jpeg. (Default: png)
        :param pixelFormat: Specifies the pixel format of the generated image. The options
        are luminance, luminanceAlpha, rgb or rgba for
        PNG output and rgb or luminance for JPEG output
        :param tonemapMethod:
        :param gamma:
        :param exposure:
        :param key:
        :param burn:
        :param banner:
        :param highQualityEdges:
        :param rfilter:
        :return:
        """
        ldrfilm_props = mitcor.Properties('ldrfilm')
        ldrfilm_props['width'] = width
        ldrfilm_props['height'] = height
        ldrfilm_props['fileFormat'] = fileFormat
        if cropOffsetX is not None:
            ldrfilm_props['cropOffsetX'] = cropOffsetX
            ldrfilm_props['cropOffsetY'] = cropOffsetY
            ldrfilm_props['cropWidth'] = cropWidth
            ldrfilm_props['cropHeight'] = cropHeight
        ldrfilm_props['pixelFormat'] = pixelFormat
        ldrfilm_props['tonemapMethod'] = tonemapMethod
        ldrfilm_props['banner'] = banner
        ldrfilm_props['gamma'] = gamma
        ldrfilm_props['exposure'] = exposure
        ldrfilm_props['key'] = key
        ldrfilm_props['burn'] = burn
        ldrfilm_props['highQualityEdges'] = highQualityEdges
        ldrfilm = plugin_mngr.createObject(ldrfilm_props)
        if rfilter is not None:
            ldrfilm.addChild(rfilter)
        ldrfilm.configure()
        self.film = ldrfilm


class NumpyFilm(object):
    def __init__(self, width=default_width, height=default_height, cropOffsetX=None, cropOffsetY=None,
                 cropWidth=None, cropHeight=None, fileFormat='numpy', pixelFormat=default_pixelFormat,
                 digits=6, variable='data', highQualityEdges=default_highQualityEdges, rfilter=None):
        """

        :param width: Width of rendered image in pixels.
        :param height: Height of rendered image in pixels.
        :param cropOffsetX: The crop parameters can optionally be provided to select a subrectangle
        of the output. In this case, Mitsuba will only render
        the requested regions. (Default: Unused). Integer.
        :param cropOffsetY: Integer.
        :param cropWidth: Integer.
        :param cropHeight: Integer.
        :param fileFormat: Specifies the desired output format;must be one of matlab,
        mathematica, or numpy. (Default: numpy). String.
        :param pixelFormat: Specifies the desired pixel format of the generated image.
        The options are luminance, luminanceAlpha, rgb, rgba,
        spectrum, and spectrumAlpha. In the latter two cases,
        the number of written channels depends on the value assigned
        to SPECTRUM_SAMPLES during compilation (see Section
        4 for details) (Default: luminance). String input.
        :param digits: Number of significant digits to be written. Default 6.
        :param variable: Name of the variable in the numpy/matlab/mathematica file
        :param highQualityEdges:
        :param rfilter: Reconstruction filter. Default is a simple box filter (BoxFilter class)
        :return:
        """
        hdrfilm_props = mitcor.Properties('mfilm')
        hdrfilm_props['width'] = width
        hdrfilm_props['height'] = height
        hdrfilm_props['fileFormat'] = fileFormat
        if cropOffsetX is not None:
            hdrfilm_props['cropOffsetX'] = cropOffsetX
            hdrfilm_props['cropOffsetY'] = cropOffsetY
            hdrfilm_props['cropWidth'] = cropWidth
            hdrfilm_props['cropHeight'] = cropHeight
        hdrfilm_props['pixelFormat'] = pixelFormat
        hdrfilm_props['digits'] = digits
        hdrfilm_props['variable'] = variable
        hdrfilm_props['highQualityEdges'] = highQualityEdges
        hdrfilm = plugin_mngr.createObject(hdrfilm_props)
        if rfilter is not None:
            hdrfilm.addChild(rfilter)
        hdrfilm.configure()
        self.film = hdrfilm

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
        z = np.sin(altitude)  # Z-axis is towards the zenith
        hyp = np.cos(altitude)
        y = hyp * np.sin(azimuth)
        x = hyp * np.cos(azimuth)
        light_direction = (-x, -y, -z)  # Light travels in opposite direction to vector towards the sun
        self.add_directional_light(direction=light_direction, irradiance=irradiance, samplingWeight=samplingWeight)

    def add_radiant_environment_map(self, filename, saz=0.0, scale=1.0, toWorld=None, gamma=1.0, cache=None,
                                    samplingWeight=1.0):
        """
        Add a radiant environment map (envmap) to a Mitsuba scene taken from a file (typically OpenEXR format).
        :param filename: Filename of the file to fetch the REM from. Supports any image format supported by Mitsuba.
        For `MORTICIA` purposes, only use OpenEXR.
        :param saz: Solar azimuth angle in degrees. Once the REM has been rotated such that the zenith is along the
        Z-axis, it is then rotated to place the solar aureole at the given azimuth. This is a positive rotation about
        the Z-axis by the solar azimuth (saz). It is assumed that the REM is provided such that the aureole
        :param scale: Scale the REM by this scalar value. Defaults to 1.0.
        :param toWorld: Transformation (only rotation matters) of the environment map. By default, the REM will be
        rotated so that the zenith is along the X-axis, the X-axis is towards the North and the Y-axis is towards the
        East. This involves a 90 degree rotation about the X-axis.
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
        if saz != 0.0:
            toWorld.rotate((0.0, 0.0, 1.0), saz)  # Rotate by solar azimuth about the new z-axis
        envmap_props = mitcor.Properties('envmap')
        envmap_props['filename'] = self.fileResolver.resolve(filename)
        envmap_props['scale'] = scale
        envmap_props['toWorld'] = toWorld.xform
        envmap_props['gamma'] = gamma
        if cache is not None:
            envmap_props['cache'] = cache
        envmap_props['samplingWeight'] = samplingWeight
        envmap = self.plugin_mngr.createObject(envmap_props)
        envmap.configure()
        self.scene.addChild('environment', envmap)


