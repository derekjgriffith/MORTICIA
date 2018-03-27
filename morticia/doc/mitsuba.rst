The Mitsuba Rendering System
============================
`Mitsuba <http://www.mitsuba-renderer.org/>`_ is an open-source, research grade, physically-based rendering system in
 the spirit of `PBRT <http://www.pbrt.org/>`_. Mitsuba meets most of the needs of ``MORTICIA``, specifically
 - Capable of spectral radiometric rendering with an arbitrary number of spectral bins (recompilation of Mitsuba is
 required when  altering the number of spectral bins).
 - Capable of using spectral radiant environment maps with any number of spectral bins in the lat/long format. These
 environment maps can be generated using the `libRadtran`/``MORTICIA`` integration and written to OpenEXR files using
  the
 `morticia.rad.librad.RadEnv` class.
 - Option to output linear and unprocessed radiometric quantities. That is, if the input quantities are absolute
 radiometric quantities, then the output radiances are physically correctly scaled.
 - Capable of making any object an emitter with a specific spectral radiance. This is required when creating
 renderings in the thermal spectrum. A current weakness in this area is that Mitsuba cannot model the effect of
 "angular emissivity", where the emissivity of a material varies with viewing angle. This can result in very "flat"
 renderings of curved objects in the thermal spectrum.
 - Large variety of physically-based surface BSDF models.
 - Highly modular source code structure, making it possible to add integrators, BSDF models, emitters etc. in a
 relatively simple manner. For example, an emitter having radiance that depends on viewing angle could be
 implemented in order to simulate angular emissivity effects.
 - Scene file format is XML, allowing for simple generation and pre-processing of scene files with a variety of
XML tools.
 - Python interface for scene generation and full control of the rendering process.
 - Mitsuba has fully-integrated and scalable parallel computing capability that can harness compute cores across a
 private network.
 - Mitsuba has a GUI (``mtsgui``) which is very useful for a number of visualisation purposes.


Spectral Rendering
------------------
For full radiometrically accurate rendering, Mitsuba must be compiled in the spectral mode with the number of spectral
bins set to 4 or more. The optimal number of bins depends on numerous factors, such as
- If Mitsuba is to be run in parallel mode across a number of network compute resources
- The complexity of the scenes to be rendered
- The number of spectral bins that are actually required for the problem at hand

Note that Mitsuba assumes that all spectral samples lie in the range of 360 nm to 830 nm (this will be called the
"visible spectrum" assumption). Certain features of Mitsuba rely on this assumption. However, for ``MORTICIA`` the need
is to cover any spectral region in the
optical domain. For example, the black body models in Mitsuba rely on the visible spectrum assumption. Any files
written using Mitsuba film models that output RGB, sRGB, XYZ color spaces also rely on the visible spectrum assumption.
The Mitsuba black body models are avoided in MORTICIA and spectra are given without wavelengths in Mitsuba scene
descriptions. ``MORTICIA`` keeps track of the actual wavelengths represented by the Mitusba spectral bins. Output of
tonemapped data in a tristimulus color space such as RGB or XYZ can still be useful for visualisation
purposes. These are essentially "false color" images when the MORTICIA scene violates the visible spectrum
assumption. The high dynamic range film (``hdrfilm``) with XYZ tristimulus outputs is the recommened method of creating
false
color images from Mitsuba. Since this will generate color information that could be "out-of-gamut" for typical display
systems, intelligent conversion to sRGB or RGB would still be required. Alternatively, simply write directly to an
RGB false color image using the ``hdrfilm`` or ``ldrfilm`` Mitsuba output plugins.

The implications of the need to recompile Mitsuba for a specific number of spectral bins are as follows:
- Renderings for scenes having more spectral bins than in the current compilation of Mitsuba will require a batch of
Mitsuba runs, with each run dealing with a subset of spectral bins. Parallel runs of Mitsuba are possible, so this is
 not necessarily a disadvantage, but scene processing overhead and memory utilisation is increased over a compilation
 of Mitsuba with more spectral bins.
- When writing REMs to OpenEXR files for Mitsuba to use as environment maps, each EXR file must have the correct
number of channels (spectral bins).
- MORTICIA requires access to an RGB installation as well as a full spectral installation of Mitsuba on the main
host platform (the user host). Since it is difficult to compile a spectral Mitsuba for Windows, it is
recommended that the main host platform be Linux, preferably Ubuntu 12.04 LTS.
- Problems have also been encountered accessing a runnning instance of mtssrv (parallel Mitsuba render service) on
a Linux platform from a Windows platform. This is another reason to use Linux exclusively.

It is possible to establish how many spectral bins a particular installation of Mitsuba has been compiled for within
the `mtsgui` GUI utility. Go to Help->About and look for the compilation flag `SPECTRUM_SAMPLES`.

Setting up Mitsuba
------------------
Refer to the Mitsuba manual for configuration and compilation. Having selected and copied a configuration file
to the file `config.py` in the root Mitsuba directory, edit this file and look in the CXXFLAGS parameter for
'-DSPECTRUM_SAMPLES=3' and change the number of spectrum samples to 4 or higher. The number of samples (bins) to use
depends on the factors listed in the previous section. If no other configuration is required, Mitsuba is then compiled
using `scons` as specified in the manual. For general multispectral work (e.g. using the Kato parametrization), consider
using 8 spectral bins. A full Kato run (all wavelengths) has 32 bins. This would require 4 runs of Mitsuba if
compiled with 8 bins. A Kato run over the VNIR spectrum is 16 bins, which would require 2 runs.

On Linux, suppose the root folder for Mitsuba is ~/Mitsuba/mitsuba and the `bash` shell is being used. Add the following
lines to the `.bashrc` file in the home directory:

# Set up Mitsuba
export MITSUBA_PYVER=2.7
source Mitsuba/mitsuba/setpath.sh

This should give access to `mtsgui` and the `mitsuba` render command from any non-login shell.


Emitters in Mitsuba
-------------------
The `sky`, `sun` and `sunsky` emitters within Mitsuba scenes are only used in the MORTICIA context for creating
presentation images. These emitter types must be avoided for quantitative work. Instead, a REM from libRadtran is used
 within an `envmap` emitter type for the diffuse component and a `directional` emitter for the direct solar component.

 In the thermal spectrum, the `directional` emitter falls away and only the `source thermal` environment map is used.

Coordinate System in Mitsuba
----------------------------
The environment maps (``envmap`` emitter) in Mitsuba is the only place in the documentation where the coordinate
system is mentioned in absolute terms. This coordinate system has +Y towards the zenith and -Y to nadir. The +X
direction is towards the left when viewing with +Y upwards. The natural (topocentric) coordinate system for ``MORTICIA``
is with +Z towards the zenith and -Z at nadir. +X is towards the north and +Y towards the east, giving a
left-handed coordinate system. In the broader context, the earth-centered, earth-fixed (ECEF) system also known as
the earth-centered rotational (ECR) coordinate system is right-handed with +Z towards the north pole, +X
through the prime meridian (Greenwich) and +Y through 90 degrees longitude measured positive east from the prime
merdian.

The Mitsuba world coordinate system is left-handed (as for PBRT) +Y typically towards the zenith, while the ``MORTICIA``
coordinate system is left-handed with +Z towards the zenith.
A coordinate transform is therefore required whe moving from ``MORTICIA`` coordinates to Mitsuba world
coordinates. The recommended method is to transform the REM coordinates in Mitsuba so that the +Z axis is upward by
rotating +90 degrees about the x-axis.
This is typically as follows::

<emitter type="envmap" >
 <string name="filename" value="REMfromMORTICIA.exr"/>
 <transform name="toWorld">
  <rotate x="1" angle="90"/>
 </transform>
</emitter>

Target models (vehicles, personnel etc.) should be edited so that they are orientated with +z upwards and such that the
normal direction of travel is +x (north in MORTICIA space). The model should be implemented as a shapegroup in the
Mitsuba scene file. This allows for orientation and placement of the entire target model in Mitsuba world coordinates.
It also allows for multiple instances of the target to be created in the Mitsuba scene file at lower computational cost.

The default origin of the world topocentric coordinate system is assumed to be a point at sea level and at nadir from
the location of the sensor. That is, the x-coordinate and y-coordinate are zero and the z-coordinate is equal to the
height (altitude) of the sensor above mean sea level (AMSL).

Sensor transformations in the Mitusba scene file should use the *lookat* form.


Mitsuba Integrators
------------------

Mitsuba has a number of "integrators", being the plugins that actually implement different rendering schemes.
The path tracer (``path``) is the integrator to be selected for general purposes, where there is direct and indirect
illumination of the scene. For faster renders, the direct illumination integrator (``direct``) provides good quality
renders without indirect illumination components. For the path tracer, the Hammersley QMC or Halton sampler is
preferred, with as many as 256 samples per pixel or more to reduce monte carlo noise.

If volumetric (participating) media such as smoke or fog are involved, the extended volumetric path tracer
(``volpath``) should be considered.

Mitsuba and Atmospheric Effects
-------------------------------
Mitsuba does not compute the effects of atmospheric scattering and absorption in the scene ("participating media") by
 default. While it is possible to define such media and use a Mitsuba integrator (e.g. volumetric integrator) that
 takes such effects into account, it would be a very difficult approach to use in this case. The model would somehow
 have to be harmonised with the atmospheric model using in the RT code (libRadtran in the case of ``MORTICIA``).

 Mitsuba is used in the most simple cases for rendering "at-target" radiance of a relatively small target (aircraft,
 vehicle or man) inserted into a radiant environment computed with libRadtran. If the scene contains objects at
 multiple distances from the sensor or if scene elements are very large (a ground plane for example) then special
 measures must be taken to model the at-sensor radiance correctly. When the sensor is at sufficient distance from the
  target element, then the effects of atmospheric path radiance and absorption must be modelled.

Ground Planes and Terrain
=========================
Rather than insert the ground plane into the Mitsuba scene, the best approach for uniform ground is simply to allow
the environment map (REM computed using libRadtran) to do the work. However, in more realistic simulations, the
ground plane will have spatial variations of reflectance and/or temperature. The first order approximation for
handling such scenes is to the use the so-called Independent Pixel Approximation (IPA). In this approach, the REMs
are computed for a range of ground reflectance and/or temperature values. REMs for any spectral surface reflectance
or temperature can then be interpolated from the REM data. Surface reflectance or temperature does not directly
influence path transmittance, so transmittance calculations are not repeated for each reflectance ("albedo" in
libRadtran). This method only deals with lambertian surfaces. Introduction of surface BRDF complicates matters greatly.

The IPA can be used for flat ground planes to perform first order approximate renderings for scenes having a
flat lambertion ground surface that is spatially non-uniform in terms of diffuse reflectance (spectral and spatial
non-uniformity).

Terrain can also be accommodated in the first-order IPA approach using a depth map to the terrain from the sensor.
Mitsuba can be used to compute this depth map, provided that it has been compiled in the RGB (non-spectral) mode.
Mitsuba does have a height map geometry shape that can be used for modelling terrain.

When compiling REMs for the solar spectrum it is recommened to compute using at least 3 albedo values (0, 0.5 and 1
.0), since path radiances are not exactly linear with surface reflectance. REMs for any surface reflectance in any
spectral bin is computed using multi-dimensional linear interpolation from this REM dataset.

In the thermal spectrum, since radiance is not not linear with temperature, it may be necessary to compute the
environment map with a significant number of surface temperatures. Alternatively, since radiance is generally
porportional to
 the fourth power of the temperature, a 4th order polynomial interpolation scheme could be used with fewer
 temperatures. However, since thermal REMs have no azimuthal dependence, it is much less costly to compute them in
 the first instance compared to solar spectrum REMs, so increasing the number of temperatures is not that costly.

General Notes on Mitsuba
------------------------

As with libRadtran, Mitsuba is not provided with ``MORTICIA``. Those wishing to use the capabilities of libRadtran or
Mitsuba will have to download, compile and install those packages on any required compute platforms and set up
supporting libraries. Correct usage of libRadtran and Mitsuba require significant insight into the relevant
knowledge domains. An effort is made to provide reasonable defaults for the many inputs that these packages require.

The Mitsuba GUI (`mtsgui`) can only read OpenEXR files with more than 3 channels if compiled with the
SPECTRAL_SAMPLES flag set higher than 3.

Mitsuba Workflows
=================
The organisation of a scene in Mitsuba follows the general practice within PBRT-like rendering systems. The Mitsuba
scene definition file is eXtensible Markup Language (XML) and the manual provides details on how to define the various
elements of a scene. One of the fastest ways to load geometry is to convert all shapes into one or more Mitsuba
serialized geometry files (extension `.serialized`). The .xml scene file then accesses any number of shapes from these
binary files. Geometry can be created using a large variety of tools or converted from `.obj` or `.dae` (Collada)
format. The Blender application can also be used to both create and convert geometry elements to either of these
formats, which can then be converted to Mitsuba format using the `mtsgui` tool provided with Mitsuba.

The very simple `.stl` format can be used for plain geometry creation and import into Blender. Many CAD applications
such as DesignSpark Mechanical can export `.stl` format files. `FreeCad <https://www.freecadweb.org/>`_ is another free
and open source tool that can be used to create geometry or convert geometry to or from `.step`, `.dae`, `.obj`, `.stl`
or other file formats.

Once a collection of indexed shapes are available in a `.serialized` file, the contents of the `.xml` scene file will
reference the shapes by index and apply BSDFs, textures, radiance and other properties to the shapes. Some properties
(e.g. texture) can be varied within a shape by using UV coordinate mapping. Simple shape geometry (spheres, cubes,
cylinders etc.) can be created within the `.xml` scene file, but complex geometry is best contained in the
native Mitsuba `.serialized` format.

Blender
-------

The `Blender <https://www.blender.org/>`_ application can be used for 3D model imports, editing, texturing and
exporting to Collada (`.dae`) or
`.obj` file formats which can then be imported into `Mitsuba`. `Blender` is a very capable environment for these
purposes, but has a complex and unique user interface together with a steep learning curve. Background knowledge
with respect to meshes, UV-mapping and texturing are generally required.

Restructuring of 3D models or renaming of model components is best done in Blender, or in the originating CAD software.

The Collada (`.dae`) exporter may still have a number of bugs. Missing or displaced components can occur, especially
 in the case of animation components. Joining components into single parts based on BSDFs can reduce the amount
 of editing required on the imported `Mitsuba` model.

`Blend Swap <http://www.blendswap.com/>`_ can be a useful source of Creative Commons 3D models.

Another useful tool for mesh visualisation, texturing, analysis, repair and format conversion is
`MeshLab <http://www.meshlab.net/>`_.

Sourcing of Models
------------------
Models can be obtained from a variety of online sources. There is large variance in the quality and organisation
of these models. Many free and commercial models are available (e.g. search `Yobi3D <http://www.yobi3d.com/>`_).
The licencing of the models must be carefully observed. Commercial models may not be distributed.
Many models are available under Creative Commons licencing, which has various levels of usage restriction.

A large repository of mainly commercial models can be found at `TurboSquid <https://www.turbosquid.com/>`_

MORTICIA Interaction with Mitsuba
---------------------------------
Mitsuba provides a powerful Python API which exposes most of the functionality of the C++ API. Documentation can be
found at the `Mitsuba Python bindings <http://www.mitsuba-renderer.org/api/group__libpython.html>`_ web page.
This is the primary way in which MORTICIA interacts with Mitsuba. An alternative method is to interact via Mitsuba
`.xml` scene files using the `lxml` package. The Jupyter notebook examples will show interaction via the Python API.

Transforming Mitsuba Scene Files
--------------------------------
Mitsuba scene files can be transformed using `XSLT` and queried using `XQuery`. The Python package `lxml` can perform
`XSLT` transformations. Microsoft Visual Studio Code is a free and powerful editor system and there is an extension
called `XML Tools` by Josh Johnson, which allows experimentation with `XQuery` and `XPath`.
The `BaseX` XML database system can also be used to manage XML documents and execute `XQuery` queries.

Viewing Mitsuba Outputs
-----------------------
In full, High Dynamic Range (HDR) spectral radiometric mode, `Mitsuba` output radiances to an OpenEXR file with
multiple (>4)
spectral channels. These files can be viewed in `mtsgui`, but only if the OpenEXR file has exactly the number of
channels for which `Mitsuba` (and therefore `mtsgui`) has been compiled. A more general OpenEXR viewer, which
allows channel selection as well gain and gamma adjustments for easier viewing, is
`mrViewer <http://mrviewer.sourceforge.net/>`_. While `IrfanView <http://mrviewer.sourceforge.net/`_ is a popular
and useful general image viewer, it cannot deal with HDR OpenEXR files from `Mitsuba`.

The recommended OpenEXR viewer for use in conjunction with `MORTICIA` and `Mitsuba` is therefore
`mrViewer <http://mrviewer.sourceforge.net/>`_. It can also be used to view Radiant Environment Maps (REM)
calculated by `libRadtran`/`MORTICIA`.



