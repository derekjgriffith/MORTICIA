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

class AnimatedTransform(object):
    """
    This class encapsulates animation transformations for Mitsuba objects (mainly target geometry components and
    sensor locations). Animations provide a series of transforms at a sequence of times. The time is assumed to be in
    seconds relative to the epoch (start time) of the scenario. The animation transformation is then effectively a
    list of Transform objects with associated elapsed time in seconds. Mitsuba does not define the time (or any other)
    physical units. These are left to the user. However, the canonical MORTICIA units for time in Mitsuba animations
    are seconds.
    """
    def __init__(self, time=0.0, xform=None):
        self.xform = mitcor.AnimatedTransform()
        if xform is not None:
            self.xform.appendTransform(time, xform.xform)

    def appendTransform(self, time, xform):
        self.xform.appendTransform(time, xform.xform)

# Classes for Mitsuba scene geometry shapes


class Shape(object):
    """
    In Mitsuba, shapes define surfaces that mark transitions between different types of materials. For
    instance, a shape could describe a boundary between air and a solid object, such as a piece of rock.
    Alternatively, a shape can mark the beginning of a space that contains a participating medium such as smoke or
    steam. Finally, a shape can be used to create an object that emits light on its own, such as a black body emitter.
    """
    shape_counter = 0

    def __init__(self, shape_props, toWorld=None,  bsdf=None, emit=None, iden=None):
        self.shape = plugin_mngr.createObject(shape_props)
        if iden is not None:
            self.iden = iden
        else:
            self.iden = 'shape_' + self.type + '_' + str(Shape.shape_counter)
            Shape.shape_counter += 1
        if toWorld is not None:
            shape_props['toWorld'] = toWorld.xform
        if bsdf is not None:
            # Should check here if bsdf is a reference to declared BSDF
            self.shape.addChild(bsdf.bsdf_type, bsdf.bsdf)
        if emit is not None:
            self.shape.addChild(emit.emit_type, emit.emitter)
        self.shape.setID(self.iden)
        self.shape.configure()

    def __str__(self):
        return str(self.shape)


class ShapeCube(Shape):
    def __init__(self, toWorld=None, flipNormals=False, bsdf=None, emit=None, iden=None):
        self.type = 'cube'
        shape_props = mitcor.Properties(self.type)
        shape_props['flipNormals'] = flipNormals
        super(ShapeCube, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeSphere(Shape):
    def __init__(self, center=(0.0, 0.0, 0.0), radius=1.0, toWorld=None, flipNormals=False, bsdf=None, emit=None,
                 iden=None):
        self.type = 'sphere'
        shape_props = mitcor.Properties(self.type)
        shape_props['center'] = mitcor.Point(center[0], center[1], center[2])
        shape_props['radius'] = radius
        shape_props['flipNormals'] = flipNormals
        super(ShapeSphere, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeCylinder(Shape):
    def __init__(self, toWorld=None, p0=(0.0, 0.0, 0.0), p1=(0.0, 0.0, 1.0), radius=1.0, flipNormals=False, bsdf=None,
                 emit=None, iden=None):
        self.type = 'cylinder'
        shape_props = mitcor.Properties(self.type)
        shape_props['p0'] = mitcor.Point(p0[0], p0[1], p0[2])
        shape_props['p1'] = mitcor.Point(p1[0], p1[1], p1[2])
        shape_props['radius'] = radius
        shape_props['flipNormals'] = flipNormals
        super(ShapeCylinder, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeRectangle(Shape):
    def __init__(self, toWorld=None, flipNormals=False, bsdf=None, emit=None, iden=None):
        self.type = 'rectangle'
        shape_props = mitcor.Properties(self.type)
        shape_props['flipNormals'] = flipNormals
        super(ShapeRectangle, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeWavefrontOBJ(Shape):
    def __init__(self, filename, toWorld=None, faceNormals=False, maxSmoothAngle=None, flipNormals=False,
                 flipTexCoords=True, collapse=False, bsdf=None, emit=None, iden=None):
        self.type = 'obj'
        shape_props = mitcor.Properties(self.type)
        shape_props['filename'] = filename
        shape_props['faceNormals'] = faceNormals
        if maxSmoothAngle is not None:
            shape_props['maxSmoothAngle'] = maxSmoothAngle
        shape_props['flipNormals'] = flipNormals
        shape_props['maxSmoothAngle'] = maxSmoothAngle
        shape_props['flipTexCoords'] = flipTexCoords
        shape_props['collapse'] = collapse
        super(ShapeWavefrontOBJ, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeStanfordTriangles(Shape):
    def __init__(self, filename, toWorld=None, faceNormals=False, maxSmoothAngle=None, flipNormals=False,
                 srgb=True, bsdf=None, emit=None, iden=None):
        self.type = 'ply'
        shape_props = mitcor.Properties(self.type)
        shape_props['filename'] = filename
        shape_props['faceNormals'] = faceNormals
        if maxSmoothAngle is not None:
            shape_props['maxSmoothAngle'] = maxSmoothAngle
        shape_props['flipNormals'] = flipNormals
        shape_props['maxSmoothAngle'] = maxSmoothAngle
        shape_props['srgb'] = srgb
        super(ShapeStanfordTriangles, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeSerialized(Shape):
    def __init__(self, filename, toWorld=None, shapeIndex=0, faceNormals=False, maxSmoothAngle=None,
                 flipNormals=False, bsdf=None, emit=None, iden=None):
        self.type = 'serialized'
        shape_props = mitcor.Properties(self.type)
        shape_props['shapeIndex'] = shapeIndex
        shape_props['filename'] = filename
        shape_props['faceNormals'] = faceNormals
        if maxSmoothAngle is not None:
            shape_props['maxSmoothAngle'] = maxSmoothAngle
        shape_props['flipNormals'] = flipNormals
        shape_props['maxSmoothAngle'] = maxSmoothAngle
        super(ShapeSerialized, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeGroup(Shape):
    def __init__(self, id, shape_list):
        self.type = 'shapegroup'
        self.id = id
        shape_props = mitcor.Properties(self.type)
        shape_props['id'] = id
        shape = plugin_mngr.createObject(shape_props)
        for shape_child in shape_list:
            shape.addChild(shape_child.shape_type + '_' + shape_child.id, shape_child.shape)
        shape.configure()
        self.shape = shape


class ShapeInstance(Shape):
    def __init__(self, shapegroup, id='noid', toWorld=None):
        self.type = 'instance'
        self.id = id
        shape_props = mitcor.Properties(self.type)
        if toWorld is not None:
            shape_props['toWorld'] = toWorld.xform
        shape = plugin_mngr.createObject(shape_props)
        shape.addChild(shapegroup)
        shape.configure()
        self.shape = shape


class ShapeHair(Shape):
    def __init__(self, filename, toWorld=None, radius=0.025, angleThreshold=1, reduction=0.0, bsdf=None, emit=None,
                 iden=None):
        self.type = 'hair'
        shape_props = mitcor.Properties(self.type)
        shape_props['filename'] = filename
        shape_props['angleThreshold'] = angleThreshold
        shape_props['reduction'] = reduction
        super(ShapeHair, self).__init__(shape_props, toWorld, bsdf, emit, iden)


class ShapeHeightField(Shape):
    def __init__(self, toWorld=None, shadingNormals=True, flipNormals=False, width=100, height=100, scale=1.0,
                 filename=None, texture=None, bsdf=None, emit=None, iden=None):
        self.type = 'heightfield'
        shape_props = mitcor.Properties(self.type)
        shape_props['shadingNormals'] = shadingNormals
        shape_props['flipNormals'] = flipNormals
        shape_props['scale'] = scale
        if texture is not None and filename is None:
            shape_props['width'] = width
            shape_props['height'] = height
            shape_props['texture'] = texture
        elif filename is not None and texture is None:
            shape_props['filename'] = filename
        else:
            warnings.warn('heightfield must be either procedural or read from a file')
        super(ShapeHeightField, self).__init__(shape_props, toWorld, bsdf, emit, iden)

# End of classes for Mitsuba geometrical shapes

# Class for Mitsuba spectra (reflectance, radiance etc.)
# Much more work needed here to tie up with MORTICIA Basis functions and color spaces

class Spectrum(object):
    """
    Encapsulates spectral vector type of data, which may provide reflectance, radiance, irradiance or perhaps other
    functions of light wavelength.
    """
    def __init__(self, spectraldata, intent='reflectance', units=None, basis=None):
        self.intent = intent
        # The spectrum must provide a two-column numpy array, where the first column is the wavelength in canonical
        # units of nm and the second column is the spectral data (reflectance, radiance, etc.)
        if isinstance(spectraldata, np.ndarray):
            self.wvl = spectraldata[:,0]
            self.data = spectraldata[:,1]
            self.spectrum = mitcor.Spectrum(self.data.tolist())
        elif isinstance(spectraldata, float):
            self.spectrum = mitcor.Spectrum(spectraldata)
            self.wvl = np.nan
            self.data = spectraldata
        elif isinstance(spectraldata, list):
            if len(spectraldata) != SPECTRUM_SAMPLES:
                warnings.warn('Provided spectrum length does not match Mitsuba SPECTRUM_SAMPLES length.')
            self.wvl = np.nan
            self.spectrum = mitcor.Spectrum(spectraldata)
            self.data = np.array(spectraldata)
        self.units = units
        self.basis = basis


# Classes for Mitsuba surface scattering models (BSDF)

class Bsdf(object):
    """
    Surface scattering models describe the manner in which light interacts with surfaces in the scene.
    They conveniently summarize the mesoscopic scattering processes that take place within the material
    and cause it to look the way it does. This represents one central component of the material system in
    Mitsuba. Another part of the renderer concerns itself with what happens in between surface interactions.

    To achieve realistic results, Mitsuba comes with a library of both general-purpose surface scattering
    models (smooth or rough glass, metal, plastic, etc.) and specializations to particular materials (woven
    cloth, masks, etc.).

    Some model plugins fit neither category and can best be described as modifiers
    that are applied on top of one or more scattering models.
    Throughout the documentation and within the scene description language, the word BSDF is used
    synonymously with the term "surface scattering model".
    """
    bsdf_counter = 0  # Global counters to assign names (id for reference purposes) to BSDFs as they are created.
    def __init__(self, bsdf_props, texture=None, iden=None):
        if iden is None:
            self.iden = 'bsdf_' + self.type + '_' + str(Bsdf.bsdf_counter)
            Bsdf.bsdf_counter += 1
        else:
            self.iden = iden
        bsdf = plugin_mngr.createObject(bsdf_props)
        if texture is not None:
            bsdf.addChild(texture.type + '_' + texture.iden, texture.texture)
        bsdf.setID(self.iden)
        bsdf.configure()
        self.bsdf = bsdf

    def __str__(self):
        return str(self.bsdf)


class BsdfDiffuse(Bsdf):
    def __init__(self, reflectance=None, texture=None, iden=None):
        self.type = 'diffuse'
        bsdf_props = mitcor.Properties(self.type)
        if reflectance is not None:
            bsdf_props['reflectance'] = reflectance
        super(BsdfDiffuse, self).__init__(bsdf_props, texture, iden)


class BsdfRoughDiffuse(Bsdf):
    def __init__(self, reflectance=None, alpha=None, useFastApprox=False, texture=None, iden=None):
        self.type = 'roughdiffuse'
        props = mitcor.Properties(self.type)
        if reflectance is not None:
            props['reflectance'] = reflectance
        if alpha is not None:
            props['alpha'] = alpha
        props['useFastApprox'] = useFastApprox
        super(BsdfRoughDiffuse, self).__init__(props, texture, iden)


class BsdfSmoothDielectric(Bsdf):
    pass


class BsdfThinDielectric(Bsdf):
    pass


class BsdfRoughDielectric(Bsdf):
    pass


class BsdfSmoothConductor(Bsdf):
    pass


class BsdfRoughCOnductor(Bsdf):
    pass


class BsdfSmoothPlastic(Bsdf):
    pass


class BsdfRoughPlastic(Bsdf):
    pass


class BsdfSmoothDielectricCoating(Bsdf):
    pass


class BsdfRoughDielectricCoating(Bsdf):
    pass


class BsdfBumpMap(Bsdf):
    pass


class BsdfAnisotropicWard(Bsdf):
    pass


class BsdfMixture(Bsdf):
    pass


class BsdfBlended(Bsdf):
    pass


class BsdfTwoSided(Bsdf):
    pass


class BsdfDiffuseTransmitter(Bsdf):
    pass


class BsdfHanrahanKruege(Bsdf):
    pass


class BsdfIrawanMarschner(Bsdf):
    pass



# End of BSDF classes

# Classes for Mitsuba textures


class Texture(object):
    """
    In Mitsuba, textures are objects
    that can be attached to certain surface scattering model parameters to introduce spatial variation.
    In the documentation, these are listed as supporting the "texture" type.
    """
    texture_counter = 0
    def __init__(self, texture_props, iden=None):
        if iden is None:
            self.iden = 'texture_' + self.type + '_' + str(Texture.texture_counter)
            Texture.texture_counter += 1
        else:
            self.iden = iden
        texture = plugin_mngr.createObject(texture_props)
        texture.setID(self.iden)
        texture.configure()
        self.texture = texture

    def __str__(self):
        return str(self.texture)


class TextureBitmap(Texture):
    def __init__(self, filename, wrapMode='repeat', wrapModeU=None, wrapModeV=None, gamma=None, filterType='ewa',
                 maxAnisotropy=20.0, cache=None, uoffset=0.0, voffset=0.0, uscale=1.0, vscale=1.0, channel=None,
                 iden=None):
        self.type='bitmap'
        props = mitcor.Properties('bitmap')
        props['filename'] = filename
        props['wrapMode'] = wrapMode
        if wrapModeU is not None:
            props['wrapModeU'] = wrapModeU
        if wrapModeV is not None:
            props['wrapModeV'] = wrapModeV
        if gamma is not None:
            props['gamma'] = gamma
        props['maxAnisotropy'] = maxAnisotropy
        props['filterType'] = filterType
        if cache is not None:
            props['cache'] = cache
        props['uoffset'] = uoffset
        props['voffset'] = voffset
        props['uscale'] = uscale
        props['vscale'] = vscale
        if channel is not None:
            props['channel'] = channel
        super(TextureBitmap, self).__init__(props, iden)


class TextureCheckerboard(Texture):
    def __init__(self, color0=Spectrum(0.4), color1=Spectrum(0.2), uoffset=0.0, voffset=0.0,
                 uscale=1.0, vscale=1.0, iden=None):
        self.type = 'checkerboard'
        props = mitcor.Properties(self.type)
        props['color0'] = color0.spectrum
        props['color1'] = color1.spectrum
        props['uoffset'] = uoffset
        props['voffset'] = voffset
        props['uscale'] = uscale
        props['vscale'] = vscale
        super(TextureCheckerboard, self).__init__(props, iden)


class TextureGrid(Texture):
    def __init__(self, color0=Spectrum(0.2), color1=Spectrum(0.4), lineWidth=0.01, uscale=1.0, vscale=1.0,
                 uoffset=0.0, voffset=0.0, iden=None):
        self.type = 'gridtexture'
        props = mitcor.Properties(self.type)
        props['color0'] = color0.spectrum
        props['color1'] = color1.spectrum
        props['lineWidth'] = lineWidth
        props['uoffset'] = uoffset
        props['voffset'] = voffset
        props['uscale'] = uscale
        props['vscale'] = vscale
        super(TextureGrid, self).__init__(props, iden)


class TextureScale(Texture):
    def __init__(self, texture_in, scale):
        self.type='scale'
        props = mitcor.Properties(self.type)
        props['scale'] = scale
        texture = plugin_mngr.createObject(props)
        texture.addChild(texture_in.texture)
        texture.configure()
        self.texture = texture


class TextureVertexColors(Texture):
    pass  # To be implemented


class TextureWireframe(Texture):
    def __init__(self, interiorColor=Spectrum(0.5), edgeColor=Spectrum(0.5), lineWidth=None, stepWidth=0.5, iden=None):
        self.type = 'wireframe'
        props = mitcor.Properties(self.type)
        props['interiorColor'] = interiorColor
        props['edgeColor'] = edgeColor
        if lineWidth is not None:
            props['lineWidth'] = lineWidth
        props[stepWidth] = stepWidth
        super(TextureWireframe, self).__init__(props, iden)


class TextureCurvature(Texture):
    def __init__(self, curvature='mean', scale=1.0, iden=None):
        self.type = 'curvature'
        props = mitcor.Properties(self.type)
        props['curvature'] = curvature
        props['scale'] = scale
        super(TextureCurvature, self).__init__(props, iden)

# End of Texture classes for Mitsuba

# Emitter classes for Mitsuba


class Emitter(object):
    """
    Mitsuba supports a wide range of emitters/light sources, which can be classified into two main categories:
    emitters which are located somewhere within the scene, and emitters that surround the scene
    to simulate a distant environment.
    """
    emitter_counter = 0

    def __init__(self, emitter_props, toWorld=None, iden=None):
        if iden is None:
            self.iden = 'emitter_' + self.type + '_' + str(Emitter.emitter_counter)
            Emitter.emitter_counter += 1
        else:
            self.iden = iden
        if toWorld is not None:
            emitter_props['toWorld'] = toWorld
        emitter = plugin_mngr.createObject(emitter_props)
        emitter.setID(self.iden)
        emitter.configure()
        self.emitter = emitter

    def __str__(self):
        return str(self.emitter)

class EmitterPoint(Emitter):
    """
    This emitter plugin implements a simple point light source, which uniformly radiates illumination
    into all directions.
    """
    def __init__(self, toWorld=None, position=(0.0, 0.0, 0.0), intensity=Spectrum(1.0), samplingWeight=1.0, iden=None):
        self.type = 'point'
        props = mitcor.Properties(self.type)
        props['position'] = mitcor.Point(position[0], position[1], position[2])
        props['intensity'] = intensity.spectrum
        props['samplingWeight'] = samplingWeight
        super(EmitterPoint, self).__init__(props, toWorld, iden)


class EmitterArea(Emitter):
    """
    This plugin implements an area light, i.e. a light source that emits diffuse illumination from the
    exterior of an arbitrary shape. Since the emission profile of an area light is completely diffuse, it has
    the same apparent brightness regardless of the observer's viewing direction. Furthermore, since it
    occupies a nonzero amount of space, an area light generally causes scene objects to cast soft shadows.

    This is currently the only emitter type that can be used to produce thermal emissions from any object in the scene
    for modelling in the thermal spectrum.
    """

    def __init__(self, radiance=Spectrum(1.0), samplingWeight=1.0, iden=None):
        self.type = 'area'
        props = mitcor.Properties(self.type)
        props['radiance'] = radiance.spectrum
        props['samplingWeight'] = samplingWeight
        super(EmitterArea, self).__init__(props, None, iden)


class EmitterSpotlight(Emitter):
    """
    This plugin provides a spot light with a linear falloff. In its local coordinate system, the spot light is
    positioned at the origin and points along the positive Z direction. It can be conveniently reoriented
    using the lookat type of transformation.
    """
    def __init__(self, toWorld=None, intensity=Spectrum(1.0), cutoffAngle=20.0, beamWidth=None, texture=None,
                  samplingWeight=1.0, iden=None):
        self.type = 'spot'
        props = mitcor.Properties(self.type)
        props['intensity'] = intensity.spectrum
        props['cutoffAngle'] = cutoffAngle
        if beamWidth is not None:
            props['beamWidth'] = beamWidth
        props['samplingWeight'] = samplingWeight
        super(EmitterSpotlight, self).__init__(props, toWorld, iden)
        if texture is not None:
            self.emitter.addChild('texture', texture.texture)
        self.emitter.configure()


class EmitterDirectional(Emitter):
    """
    This emitter plugin implements a distant directional source, which radiates a specified power per
    unit area along a fixed direction. By default, the emitter radiates in the direction of the negative Z axis,
    which in the MORTICIA frame is toward nadir (sun in the zenith).
    This is emitter type used by MORTICIA for direct sunlight (or moonlight)
    """
    def __init__(self, toWorld=None, direction=(0.0, 0.0, -1.0), irradiance=Spectrum(1.0), samplingWeight=1.0,
                 iden=None):
        self.type = 'directional'
        props = mitcor.Properties(self.type)
        if direction is not None and toWorld is None:
            props['direction'] = mitcor.Vector(direction[0], direction[1], direction[2])
        elif toWorld is not None and direction is not None:
            warnings.warn('toWorld and direction inputs to directional emitter both provided. direction ignored.')
        props['samplingWeight'] = samplingWeight
        super(EmitterDirectional, self).__init__(props, toWorld, iden)



class EmitterCollimated(Emitter):
    """
    This emitter plugin implements a collimated beam source, which radiates a specified amount of
    power along a fixed ray. It can be thought of as the limit of a spot light as its field of view tends to
    zero.
    Such a emitter is useful for conducting virtual experiments and testing the renderer for correctness.
    By default, the emitter is located at the origin and radiates into the positive Z direction (0, 0, 1).
    This can be changed by providing a custom toWorld transformation.
    """
    def __init__(self, toWorld=None, power=Spectrum(1.0), samplingWeight=1.0, iden=None):
        self.type = 'collimated'
        props = mitcor.Properties(self.type)
        props['power'] = power.spectrum
        props['samplingWeight'] = samplingWeight
        super(EmitterCollimated, self).__init__(props, toWorld, iden)


class EmitterSky(Emitter):
    """
    This plugin provides the physically-based skylight model by Hosek and Wilkie. It can be used
    to create predictive daylight renderings of scenes under clear skies, which is useful for architectural
    and computer vision applications. The implementation in Mitsuba is based on code that was generously
    provided by the authors.
    The model has two main parameters: the turbidity of the atmosphere and the position of the sun.
    The position of the sun in turn depends on a number of secondary parameters, including the
    latitude, longitude, and timezone at the location of the observer, as well as the current year,
    month, day, hour, minute, and second.
    This model is used in MORTICIA only for fast visualisations of scenes under natural illumination and not for
    quantitative purposes. Rather, libRadtran REMs are used in conjunction with the directional and envmap emitters.
    """
    def __init__(self, toWorld=None, turbidity=3.0, albedo=Spectrum(0.15), year=2010, month=7, day=10, hour=15,
                 minute=0, second=0, latitude=35.6894, longitude=139.6917, timezone=9, sunDirection=None,
                 stretch=1.0, resolution=512, scale=1.0, samplingWeight=1.0, iden=None):
        self.type = 'sky'
        props = mitcor.Properties(self.type)
        props['turbidity'] = turbidity
        props['albedo'] = albedo
        props['year'] = year
        props['month'] = month
        props['year'] = year
        props['day'] = day
        props['hour'] = hour
        props['minute'] = minute
        props['second'] = second
        props['latitude'] = latitude
        props['longitude'] = longitude
        props['timezone'] = timezone
        if sunDirection is not None:
            props['sunDirection'] = mitcor.Vector(sunDirection[0], sunDirection[1], sunDirection[2])
        props['stretch'] = stretch
        props['resolution'] = resolution
        props['scale'] = scale
        props['samplingWeight'] = samplingWeight
        super(EmitterSky, self).__init__(props, toWorld, iden)

class EmitterSun(Emitter):
    """
    This plugin implements the physically-based sun model proposed by Preetham et al. [38]. Using
    the provided position and time information (see sky for details), it can determine the position of the
    sun as seen from the position of the observer.The radiance arriving at the earth surface is then found
    based on the spectral emission profile of the sun and the extinction cross-section of the atmosphere
    (which depends on the turbidity and the zenith angle of the sun).
    """
    def __init__(self, turbidity=3.0, year=2010, month=7, day=10, hour=15,
                 minute=0, second=0, latitude=35.6894, longitude=139.6917, timezone=9, sunDirection=None,
                 resolution=512, scale=1.0, sunRadiusScale=None, samplingWeight=1.0, iden=None):
        self.type = 'sun'
        props = mitcor.Properties(self.type)
        props['turbidity'] = turbidity
        props['year'] = year
        props['month'] = month
        props['year'] = year
        props['day'] = day
        props['hour'] = hour
        props['minute'] = minute
        props['second'] = second
        props['latitude'] = latitude
        props['longitude'] = longitude
        props['timezone'] = timezone
        if sunDirection is not None:
            props['sunDirection'] = mitcor.Vector(sunDirection[0], sunDirection[1], sunDirection[2])
        props['resolution'] = resolution
        props['scale'] = scale
        if sunRadiusScale is not None:
            props['sunRadiusScale'] = sunRadiusScale
        props['samplingWeight'] = samplingWeight
        super(EmitterSun, self).__init__(props, None, iden)


class EmitterSunSky(Emitter):
    """
    This convenience plugin has the sole purpose of instantiating `sun` and `sky` and merging them into
    a joint environment map. Please refer to these plugins individually for more details. This is not really
    useful in MORTICIA, since it is not possible to rotate the sky component into the correct orientation. Rather
    use `sun` and `sky` individually.
    """
    def __init__(self, turbidity=3.0, albedo=Spectrum(0.15), year=2010, month=7, day=10, hour=15,
                 minute=0, second=0, latitude=35.6894, longitude=139.6917, timezone=9, sunDirection=None,
                 stretch=1.0, resolution=512, sunScale=1.0, skyScale=1.0, sunRadiusScale=None, samplingWeight=1.0, \
                                                                                                         iden=None):
        self.type = 'sky'
        props = mitcor.Properties(self.type)
        props['turbidity'] = turbidity
        props['albedo'] = albedo
        props['year'] = year
        props['month'] = month
        props['year'] = year
        props['day'] = day
        props['hour'] = hour
        props['minute'] = minute
        props['second'] = second
        props['latitude'] = latitude
        props['longitude'] = longitude
        props['timezone'] = timezone
        if sunDirection is not None:
            props['sunDirection'] = mitcor.Vector(sunDirection[0], sunDirection[1], sunDirection[2])
        props['stretch'] = stretch
        props['resolution'] = resolution
        props['sunScale'] = sunScale
        props['skyScale'] = skyScale
        if sunRadiusScale is not None:
            props['sunRadiusScale'] = sunRadiusScale
        props['samplingWeight'] = samplingWeight
        super(EmitterSunSky, self).__init__(props, None, iden)


class EmitterEnvironmentMap(Emitter):
    """
    This plugin provides a HDRI (high dynamic range imaging) environment map, which is a type of
    light source that is well-suited for representing natural illumination. Together with the `directional` emitter
    this is the most important emitter class for MORTICIA purposes and the the REMs generated using libRadtran are used
    as the OpenEXR input.
    """
    def __init__(self, filename, toWorld=None, scale=1.0, gamma=1.0, cache=None, samplingWeight=1.0, iden=None):
        self.type = 'envmap'
        props = mitcor.Properties(self.type)
        props['filename'] = filename
        props['scale'] = scale
        props['gamma'] = gamma
        if cache is not None:
            props['cache'] = cache
        props['samplingWeight'] = samplingWeight
        super(EmitterEnvironmentMap, self).__init__(props, toWorld, iden)



class EmitterConstantEnvironment(Emitter):
    """
    This plugin implements a constant environment emitter, which surrounds the scene and radiates
    diffuse illumination towards it. This is often a good default light source when the goal is to visualize
    some loaded geometry that uses basic (e.g. diffuse) materials.
    """
    def __init__(self, radiance=Spectrum(1.0), samplingWeight=1.0, iden=None):
        self.type = 'constant'
        props = mitcor.Properties(self.type)
        props['radiance'] = radiance.spectrum
        props['samplingWeight'] = samplingWeight
        super(EmitterConstantEnvironment, self).__init__(props, None, iden)

# End of classes for Mitsuba emitters

# Classes for Mitsuba Sensors

class Sensor(object):
    """
    In Mitsuba, sensors, along with a film, are responsible for recording radiance measurements in some
    usable format. This includes default choices such as perspective or orthographic cameras, as well as
    more specialized sensors that measure the radiance into a given direction or the irradiance received
    by a certain surface.
    """
    sensor_counter = 0

    def __init__(self, sensor_props, toWorld=None, iden=None):
        if iden is None:
            self.iden = 'sensor_' + self.type + '_' + str(Sensor.sensor_counter)
            Sensor.sensor_counter += 1
        else:
            self.iden = iden
        if toWorld is not None:
            sensor_props['toWorld'] = toWorld
        sensor = plugin_mngr.createObject(sensor_props)
        sensor.setID(self.iden)
        sensor.configure()
        self.sensor = sensor

    def __str__(self):
        return str(self.sensor)


class SensorPerspective(Sensor):
    """
    This plugin implements a simple idealized perspective camera model, which has an infinitely small
    aperture. This creates an infinite depth of field, i.e. no optical blurring occurs. The camera can be
    animated to move during an exposure, hence temporal blur is still possible.
    By default, the camera's field of view is specified using a 35 mm film equivalent focal length, which is
    first converted into a diagonal field of view and subsequently applied to the camera. This assumes that
    the film's aspect ratio matches that of 35mm film (1.5:1), though the parameter still behaves intuitively
    when this is not the case. Alternatively, it is also possible to specify a field of view in degrees along a
    given axis (see the fov and fovAxis parameters).
    The exact camera position and orientation is most easily expressed using the lookAt form of transformation.
    In MORTICIA, the field of view of the camera is specified using the `fov` (in degrees) and `fovAxis` parameters.
    The default `farClip` plane is 100 000 Mitsuba length units, which in the MORTICIA context is 100 000 m. The near
    clip is 0.01 units, which is 1 cm in canonical MORTICIA units of metres.
    """
    def __init__(self, toWorld=None, focalLength='50mm', fov=None, fovAxis='x', shutterOpen=None, shutterClose=None,
                 nearClip=0.01, farClip=100000.0, iden=None):
        self.type = 'perspective'
        props = mitcor.Properties(self.type)
        if fov is not None:
            props['fov'] = fov
            props['fovAxis'] = fovAxis
        else:
            props['focalLength'] = focalLength
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        props['nearClip'] = nearClip
        props['farClip'] = farClip
        super(SensorPerspective, self).__init__(props, toWorld, iden)


class SensorPerspectiveThinLens(Sensor):
    """
    This plugin implements a simple perspective camera model with a thin lens at its circular aperture.
    It is very similar to the perspective plugin except that the extra lens element permits rendering
    with a specifiable (i.e. non-infinite) depth of field. To configure this, it has two extra parameters
    named `apertureRadius` and `focusDistance`, which are mandatory inputs.
    """

    def __init__(self, apertureRadius, focusDistance, toWorld=None,
                 focalLength='50mm', fov=None, fovAxis='x', shutterOpen=None, shutterClose=None,
                 nearClip=0.01, farClip=100000.0, iden=None):
        self.type = 'thinlens'
        props = mitcor.Properties(self.type)
        props['apertureRadius'] = apertureRadius
        props['focusDistance'] = focusDistance
        if fov is not None:
            props['fov'] = fov
            props['fovAxis'] = fovAxis
        else:
            props['focalLength'] = focalLength
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        props['nearClip'] = nearClip
        props['farClip'] = farClip
        super(SensorPerspectiveThinLens, self).__init__(props, toWorld, iden)


class SensorOrthographic(Sensor):
    """
    This plugin implements a simple orthographic camera, i.e. a sensor based on an orthographic
    projection without any form of perspective. It can be thought of as a planar sensor that measures the
    radiance along its normal direction. By default this plane is a 2 by 2 length unit (metres in MORTICIA) region in the
    XY plane (horizontal plane in MORTICIA) viewing along the +Z direction (vertically upward in MORTICIA). A
    transform can be used to scale and rotate the orthographic plane.
    """
    def __init__(self, toWorld=None, shutterOpen=None, shutterClose=None,
                 nearClip=0.01, farClip=100000.0, iden=None):
        self.type = 'orthographic'
        props = mitcor.Properties(self.type)
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        props['nearClip'] = nearClip
        props['farClip'] = farClip
        super(SensorOrthographic, self).__init__(props, toWorld, iden)


class SensorTelecentric(Sensor):
    """
    This plugin implements a simple model of a camera with a telecentric lens. This is a type of lens
    that produces an in-focus orthographic view on a plane at some distance from the sensor. Points
    away from this plane are out of focus and project onto a circle of confusion. In comparison to idealized
    orthographic cameras, telecentric lens cameras exist in the real world and find use in some
    computer vision applications where perspective effects cause problems. This sensor relates to the
    `orthographic` plugin in the same way that `thinlens` does to `perspective`.
    The configuration is identical to the orthographic plugin, except that the additional parameters
    `apertureRadius` and `focusDistance` must be provided (floats).
    """
    def __init__(self, apertureRadius, focusDistance, toWorld=None,
                 shutterOpen=None, shutterClose=None,
                 nearClip=0.01, farClip=100000.0, iden=None):
        self.type = 'telecentric'
        props = mitcor.Properties(self.type)
        props['apertureRadius'] = apertureRadius
        props['focusDistance'] = focusDistance
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        props['nearClip'] = nearClip
        props['farClip'] = farClip
        super(SensorTelecentric, self).__init__(props, toWorld, iden)


class SensorSpherical(Sensor):
    """
    The spherical camera captures the illumination arriving from all directions and turns it into a
    latitude-longitude environment map. It is best used with a high dynamic range film that has 2:1 aspect
    ratio, and the resulting output can then be turned into a distant light source using the envmap plugin.
    By default, the camera is located at the origin, which can be changed by providing a custom toWorld
    transformation.
    """
    def __init__(self, toWorld=None, shutterOpen=None, shutterClose=None, iden=None):
        self.type = 'spherical'
        props = mitcor.Properties(self.type)
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        super(SensorSpherical, self).__init__(props, toWorld, iden)


class SensorIrradianceMeter(Sensor):
    """
    This sensor plugin implements a simple irradiance meter, which measures the average incident
    power per unit area over a provided surface. Such a sensor is useful for conducting virtual experiments
    and testing the renderer for correctness. The result is normalized so that an irradiance sensor
    inside an integrating sphere with constant radiance 1 records an irradiance value of pi.

    To create an irradiance meter, instantiate the desired measurement shape and specify the sensor
    as its child. Note that when the sensor's film resolution is larger than 1 by 1, pixel will record the
    average irradiance over a rectangular part of the shape's UV parameterization.
    """
    def __init__(self, shutterOpen=None, shutterClose=None, iden=None):
        self.type = 'irradiancemeter'
        props = mitcor.Properties(self.type)
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        toWorld = None
        super(SensorIrradianceMeter, self).__init__(props, toWorld, iden)



class SensorRadianceMeter(Sensor):
    """
    This sensor plugin implements a simple radiance meter, which measures the incident power per
    unit area per unit solid angle along a certain ray. It can be thought of as the limit of a standard
    perspective camera as its field of view tends to zero. Hence, when this sensor is given a film with
    multiple pixels, all of them will record the same value.
    Such a sensor is useful for conducting virtual experiments and testing the renderer for correctness.
    By default, the sensor is located at the origin and performs a measurement in the positive Z direction
    (0, 0, 1).This can be changed by providing a custom toWorld transformation.
    """
    def __init__(self, toWorld=None, shutterOpen=None, shutterClose=None, iden=None):
        self.type = 'radiancemeter'
        props = mitcor.Properties(self.type)
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        super(SensorRadianceMeter, self).__init__(props, toWorld, iden)


class SensorFluenceMeter(Sensor):
    """
    This sensor plugin implements a simple fluence meter, which measures the average radiance passing
    through a specified position. By default, the sensor is located at the origin.
    Such a sensor is useful for conducting virtual experiments and testing the renderer for correctness.
    """
    def __init__(self, toWorld=None, shutterOpen=None, shutterClose=None, iden=None):
        self.type = 'fluencemeter'
        props = mitcor.Properties(self.type)
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        super(SensorFluenceMeter, self).__init__(props, toWorld, iden)

class SensorPerspectivePinholeDistortion(Sensor):
    """
    This plugin extends the perspective camera with support for radial distortion. It accepts an additional
    parameter named kc, which specifies the second and fourth-order terms in a polynomial
    model that accounts for pincushion and barrel distortion. This is useful when trying to match rendered
    images to photographs created by a camera whose distortion is known. When kc=0, 0, the
    model turns into a standard pinhole camera. The reason for creating a separate plugin for this feature
    is that distortion involves extra overheads per ray that users may not be willing to pay for if
    their scene doesn't use it. The `MATLAB Camera Calibration Toolbox <http:
    //www.vision.caltech.edu/bouguetj/calib_doc/>`_ by Jean-Yves Bouguet can be used to
    obtain a distortion model,
    and the first entries of the kc variable generated by this tool can directly be passed into this plugin.
    """
    def __init__(self, kc, toWorld=None, focalLength='50mm', fov=None, fovAxis='x', shutterOpen=None, shutterClose=None,
                 nearClip=0.01, farClip=100000.0, iden=None):
        self.type = 'perspective_rdist'
        props = mitcor.Properties(self.type)
        if fov is not None:
            props['fov'] = fov
            props['fovAxis'] = fovAxis
        else:
            props['focalLength'] = focalLength
        if shutterOpen is not None:
            props['shutterOpen'] = shutterOpen
        if shutterClose is not None:
            props['shutterClose'] = shutterClose
        props['nearClip'] = nearClip
        props['farClip'] = farClip
        super(SensorPerspectivePinholeDistortion, self).__init__(props, toWorld, iden)


# End of classes for Mitsuba Sensors


# Classes for different Mitsuba reconstruction filters

class Filter(object):
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
    def __init__(self, rfilter_props):
        self.rfilter = plugin_mngr.createObject(rfilter_props)
        self.rfilter.configure()
        self.radius = self.rfilter.getRadius()
        self.borderSize = self.rfilter.getBorderSize()

    def __str__(self):
        return str(self.rfilter)

class FilterBox(Filter):
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
        super(FilterBox, self).__init__(rfilter_props)


class FilterTent(Filter):
    """
    Simple tent, or triangle filter.This reconstruction filter never suffers from ringing
    and usually causes less aliasing than a naive box filter.When rendering scenes with sharp brightness
    discontinuities, this may be useful; otherwise, negative-lobed filters will be preferable (e.g.
    Mitchell-Netravali or Lanczos Sinc). The filter radius defaults to 1.0 pixels.
    """
    def __init__(self):
        self.filtertype = 'tent'
        rfilter_props = mitcor.Properties(self.filtertype)
        super(FilterTent, self).__init__(rfilter_props)


class FilterGaussian(Filter):
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
        super(FilterGaussian, self).__init__(rfilter_props)


class FilterMitchell(Filter):
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
        super(FilterMitchell, self).__init__(rfilter_props)


class FilterCatmullRom(Filter):
    """
    This is a special version of the Mitchell-Netravali filter that has
    the constants B and C adjusted to produce higher sharpness at the cost of increased susceptibility
    to ringing.
    """
    def __init__(self):
        self.filtertype = 'catmullrom'
        rfilter_props = mitcor.Properties(self.filtertype)
        super(FilterCatmullRom, self).__init__(rfilter_props)


class FilterLanczos(Filter):
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
        super(FilterLanczos, self).__init__(rfilter_props)

# End of Mitsuba reconstruction filter classes

# Classes for various times of Mitsuba sensor "film".

class Film(object):
    """
    The Film classes encapsulates Mitsuba films. From the Mitsuba manual:
    A film defines how conducted measurements are stored and converted into the final output file that
    is written to disk at the end of the rendering process. Mitsuba comes with a few films that can write
    to high and low dynamic range image formats (OpenEXR, JPEG or PNG), as well more scientifically
    oriented data formats (e.g. Numpy, MATLAB or Mathematica).
    """
    def __init__(self, film_props, rfilter=None):
        film = plugin_mngr.createObject(film_props)
        if rfilter is not None:
            film.addChild(rfilter)
        film.configure()
        self.film = film

    def __str__(self):
        return str(self.film)


class FilmHdr(Film):
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
        super(FilmHdr, self).__init__(hdrfilm_props, rfilter)

class FilmTiledHdr(Film):
    """
    This plugin implements a camera film that stores the captured image as a tiled high dynamic-range
    OpenEXR file. It is very similar to hdrfilm, the main difference being that it does not keep the
    rendered image in memory. Instead, image tiles are directly written to disk as they are being rendered,
    which enables renderings of extremely large output images that would otherwise not fit into memory
    (e.g. 100K by 100K).
    When the image can fit into memory, usage of this plugin is discouraged: due to the extra overhead
    of tracking image tiles, the rendering process will be slower, and the output files also generally do not
    compress as well as those produced by hdrfilm.
    """
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
        super(FilmTiledHdr, self).__init__(hdrfilm_props, rfilter)

class FilmLdr(Film):
    """
    This plugin implements a low dynamic range film that can write out 8-bit PNG and JPEG images
    in various configurations. It provides basic tonemapping techniques to map recorded radiance values
    into a reasonable displayable range. An alpha (opacity) channel can be written if desired. By default,
    the plugin writes gamma-corrected PNG files using the sRGB color space and no alpha channel.
    This film is a good choice when lowdynamic range output is desired and the rendering setup can be
    configured to capture the relevant portion of the dynamic range reliably enough so that the original
    HDR data can safely be discarded. When this is not the case, it may be easier to use hdrfilm along
    with the batch tonemapper.
    """
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
        super(FilmLdr, self).__init__(ldrfilm_props, rfilter)


class FilmNumpy(Film):
    """
    This plugin provides a camera film that exports spectrum, RGB, XYZ, or luminance values as a
    matrix to aMATLAB orMathematica ASCII file or a NumPy binary file.This is useful when running
    Mitsuba as simulation step as part of a larger virtual experiment. It can also come in handy when
    verifying parts of the renderer using an automated test suite.
    """
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
        super(FilmNumpy, self).__init__(hdrfilm_props, rfilter)

# Classes for Mitsuba samplers, including 3 quasi-Monte Carlo (QMC) samplers
# The Halton sampler is preferred for MORTICIA work and is the default


class Sampler(object):
    """
    When rendering an image, Mitsuba has to solve a high-dimensional integration problem that involves
    the geometry, materials, lights, and sensors that make up the scene. Because of the mathematical
    complexity of these integrals, it is generally impossible to solve them analytically. Instead, they
    solved numerically by evaluating the function to be integrated at a large number of different positions
    referred to as samples. Sample generators are an essential ingredient to this process: they produce
    points in a (hypothetical) infinite dimensional hypercube that constitute the canonical representation
    of these samples. See the Mitsuba manual for further details.
    """
    def __init__(self, sampler_props):
        sampler = plugin_mngr.createObject(sampler_props)
        sampler.configure()
        self.sampler = sampler

    def __str__(self):
        return str(self.sampler)


class SamplerIndependent(Sampler):
    """
    The independent sampler produces a stream of independent and uniformly distributed pseudorandom
    numbers. Internally, it relies on a fast SIMD version of the Mersenne Twister random number
    generator.
    This is the most basic sample generator; because no precautions are taken to avoid sample clumping,
    images produced using this plugin will usually take longer to converge. In theory, this sampler is
    initialized using a deterministic procedure, which means that subsequent runs of Mitsuba should
    create the same image. In practice, when rendering with multiple threads and/or machines, this is
    not true anymore, since the ordering of samples is influenced by the operating system scheduler.
    """
    def __init__(self, sampleCount=4):
        sampler_props = mitcor.Properties('independent')
        sampler_props['sampleCount'] = sampleCount
        super(SamplerIndependent, self).__init__(sampler_props)


class SamplerStratified(Sampler):
    """
    The stratified sample generator divides the domain into a discrete number of strata and produces a
    sample within each one of them.This generally leads to less sample clumping when compared to the
    independent sampler, as well as better convergence. Due to internal storage costs, stratified samples
    are only provided up to a certain dimension, after which independent sampling takes over.
    """
    def __init__(self, sampleCount=4, dimension=4):
        sampler_props = mitcor.Properties('stratified')
        sampler_props['sampleCount'] = sampleCount
        sampler_props['dimension'] = dimension
        super(SamplerStratified, self).__init__(sampler_props)


class SamplerLowDiscrepancy(Sampler):
    """
    This plugin implements a simple hybrid sampler that combines aspects of a Quasi-Monte Carlo sequence
    with a pseudorandom number generator based on a technique proposed by Kollig and Keller.
    It is a good and fast general-purpose sample generator and therefore chosen as the default option
    in Mitsuba. Some of the QMC samplers in the following pages can generate even better distributed
    samples, but this comes at a higher cost in terms of performance.
    """
    def __init__(self, sampleCount=4, dimension=4):
        sampler_props = mitcor.Properties('ldsampler')
        sampler_props['sampleCount'] = sampleCount
        sampler_props['dimension'] = dimension
        super(SamplerLowDiscrepancy, self).__init__(sampler_props)


class SamplerHalton(Sampler):
    """
    This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the Halton sequence.
    QMC number sequences are designed to reduce sample clumping across integration dimensions,
    which can lead to a higher order of convergence in renderings. Because of the deterministic
    character of the samples, errors will manifest as grid or moire patterns rather than randomnoise, but
    these diminish as the number of samples is increased.
    """
    def __init__(self, sampleCount=4, scramble=-1):
        sampler_props = mitcor.Properties('halton')
        sampler_props['sampleCount'] = sampleCount
        sampler_props['scramble'] = scramble
        super(SamplerHalton, self).__init__(sampler_props)


class SamplerHammersley(Sampler):
    """
    This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the Hammersley
    sequence. QMC number sequences are designed to reduce sample clumping across integration
    dimensions, which can lead to a higher order of convergence in renderings. Because of the deterministic
    character of the samples, errors will manifest as grid or moire patterns rather than randomnoise,
    but these diminish as the number of samples is increased.
    """
    def __init__(self, sampleCount=4, scramble=-1):
        sampler_props = mitcor.Properties('halton')
        sampler_props['sampleCount'] = sampleCount
        sampler_props['scramble'] = scramble
        super(SamplerHammersley, self).__init__(sampler_props)


class SamplerSobol(Sampler):
    """
    This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the Sobol sequence.
    QMC number sequences are designed to reduce sample clumping across integration dimensions,
    which can lead to a higher order of convergence in renderings. Because of the deterministic
    character of the samples, errors will manifest as grid or moire patterns rather than random noise, but
    these diminish as the number of samples is increased.
    """
    def __init__(self, sampleCount=4, scramble=-1):
        sampler_props = mitcor.Properties('halton')
        sampler_props['sampleCount'] = sampleCount
        sampler_props['scramble'] = scramble
        super(SamplerSobol, self).__init__(sampler_props)

# End of classes for Mitsuba samplers

# Classes for Mitsuba integrators



class Integrator(object):
    """
    In Mitsuba, the different rendering techniques are collectively referred to as integrators, since they
    perform integration over a high-dimensional space. Each integrator represents a specific approach
    for solving the light transport equation - usually favored in certain scenarios, but at the same time affected
    by its own set of intrinsic limitations. Therefore, it is important to carefully select an integrator
    based on user-specified accuracy requirements and properties of the scene to be rendered. See the Mitsuba manual
    for further details. In general the basic path tracer(`path`) is used for MORTICIA work unless there are specific
    reasons to use something else (e.g. scene contains participating media such as smoke).

    For faster rendering, when direct illumination is all that is required, the `direct` integrator can be used.

    The most important integrators for MORTICIA purposes are PathTracer, DirectIllumination,
    Multichannel and FieldExtraction. Other integrators will be accommodated as required.

    The AmbientOcclusion integrator can be useful for simple examination of scenes and is therefore also accommodated.
    """
    def __init__(self, integrator_props):
        integrator = plugin_mngr.createObject(integrator_props)
        integrator.configure()
        self.integrator = integrator

    def __str__(self):
        return str(self.integrator)


class IntegratorAmbientOcclusion(Integrator):
    def __init__(self, shadingSamples=1, rayLength=-1):
        self.integrator_type = 'ao'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['shadingSamples'] = shadingSamples
        integrator_props['rayLength'] = rayLength
        super(IntegratorAmbientOcclusion, self).__init__(integrator_props)


class IntegratorDirectIllumination(Integrator):
    def __init__(self, shadingSamples=1, emitterSamples=1, bsdfSamples=1, strictNormals=False, hideEmitters=False):
        self.integrator_type = 'direct'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['shadingSamples'] = shadingSamples
        integrator_props['emitterSamples'] = emitterSamples
        integrator_props['bsdfSamples'] = bsdfSamples
        integrator_props['strictNormals'] = strictNormals
        integrator_props['hideEmitters'] = hideEmitters
        super(IntegratorDirectIllumination, self).__init__(integrator_props)


class IntegratorPathTracer(Integrator):
    def __init__(self, maxDepth=-1, rrDepth=5, strictNormals=False, hideEmitters=False):
        self.integrator_type = 'path'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['maxDepth'] = maxDepth
        integrator_props['rrDepth'] = rrDepth
        integrator_props['strictNormals'] = strictNormals
        integrator_props['hideEmitters'] = hideEmitters
        super(IntegratorPathTracer, self).__init__(integrator_props)


class IntegratorSimpleVolumetric(Integrator):
    def __init__(self, maxDepth=-1, rrDepth=5, strictNormals=False, hideEmitters=False):
        self.integrator_type = 'volpath_simple'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['maxDepth'] = maxDepth
        integrator_props['rrDepth'] = rrDepth
        integrator_props['strictNormals'] = strictNormals
        integrator_props['hideEmitters'] = hideEmitters
        super(IntegratorSimpleVolumetric, self).__init__(integrator_props)


class IntegratorExtendedVolumetricPathTracer(Integrator):
    def __init__(self, maxDepth=-1, rrDepth=5, strictNormals=False, hideEmitters=False):
        self.integrator_type = 'volpath'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['maxDepth'] = maxDepth
        integrator_props['rrDepth'] = rrDepth
        integrator_props['strictNormals'] = strictNormals
        integrator_props['hideEmitters'] = hideEmitters
        super(IntegratorExtendedVolumetricPathTracer, self).__init__(integrator_props)


class IntegratorBidirectionalPathTracer(Integrator):
    def __init__(self, maxDepth=-1, rrDepth=5, lightImage=True, sampleDirect=True):
        self.integrator_type = 'bdpt'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['maxDepth'] = maxDepth
        integrator_props['rrDepth'] = rrDepth
        integrator_props['lightImage'] = lightImage
        integrator_props['sampleDirect'] = sampleDirect
        super(IntegratorBidirectionalPathTracer, self).__init__(integrator_props)


class IntegratorPhotonMap(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorProgressivePhotonMap(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorStochasticProgressivePhotonMap(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorPrimarySampleSpaceMetropolis(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorPathSpaceMetropolis(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorEnergyRedistributionPathTracer(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')

class IntegratorAdjointParticleTracer(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorAdaptive(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorVirtualPointLight(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorIrradianceCaching(Integrator):
    def __init__(self):
        warnings.warn('Not implemented yet in mortsuba.py')


class IntegratorMultichannel(Integrator):
    def __init__(self, integrator_list):
        """
        The Multichannel Integrator is a container for multiple integrators. First instantiate the sub-integrators and
        then pass them as a list to Multichannel Integrator. See the Mitsuba manual for details and examples.
        :param integrator_list: List of integrators to combine into a multichannel integrator
        Multichannel Integrators can likely not be nested i.e. do not pass a multichannel integrator as an element of
        the integrator_list.
        :return:
        """
        # Instantiate the plugin
        self.integrator_type = 'multichannel'
        integrator = plugin_mngr.create({'type': self.integrator_type})
        # Add the list of sub-integrators as children
        for integrator_child in integrator_list:
            if integrator_child.integrator_type == 'multichannel':
                warnings.warn('multichannel integrators can likely not be nested.')
            integrator.addChild(integrator_child.integrator_type, integrator_child.integrator)
        integrator.configure()
        self.integrator = integrator


class IntegratorFieldExtraction(Integrator):
    def __init__(self, field, undefined=None):
        """
        This integrator extracts a requested field from the ray intersection records of shading points and converts
        the resulting data into color values. It is meant to be used in conjunction with MultichannelIntegrator
        to dump auxiliary information (such as depth or surface normals of surfaces seen by the camera)
        into extra channels of a rendered image, for instance to create benchmark data for computer vision
        applications. Refer to the documentation of the MultichannelIntegrator for examples.
        :param field: Denotes the data that should be extracted for writing to the image data instead of radiometric
        data. The following choices are available:
            position: 3D position in world space
            relPosition: 3D position in camera space
            distance: Ray distance to the shading point
            geoNormal: Geometric surface normal
            shNormal: Shading surface normal
            uv: UV coordinate value
            albedo: Albedo value of the BSDF
            shapeIndex: Integer index of the high-level shape
            primIndex: Integer shape primitive index
        :param undefined: Value returned when there is no geometry intersection. Default is zero.
        :return:
        """
        self.integrator_type = 'field'
        integrator_props = mitcor.Properties(self.integrator_type)
        integrator_props['field'] = field
        if undefined is not None:
            integrator_props['undefined'] = undefined
        super(IntegratorFieldExtraction, self).__init__(integrator_props)

# End of classes for Mitsuba integrators

class Scene(object):
    """ The Scene class encapsulates the information that can be contained in a Mitsuba scene file.
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

    def __str__(self):
        return str(self.scene)

    def hasEnvironmentEmitter(self):
        return self.scene.hasEnvironmentEmitter()

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
        dir_light = plugin_mngr.createObject(dir_light_props)
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
        if self.hasEnvironmentEmitter():
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
        envmap = plugin_mngr.createObject(envmap_props)
        envmap.configure()
        self.scene.addChild('environment', envmap)


