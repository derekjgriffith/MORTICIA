The Mitsuba Rendering System
============================
`Mitsuba <http://www.mitsuba-renderer.org/>`_ is an open-source, research grade, physically-based rendering system in
 the spirit of `PBRT <http://www.pbrt.org/>`_. Mitsuba meets most of the needs of ``MORTICIA``, specifically
 - Capable of spectral rendering with an arbitrary number of spectral bins.
 - Capable of using spectral radiant environment maps with any number of spectral bins in the lat/long format.


Spectral Rendering
------------------
Mitsuba must be compiled in the spectral mode with the number of spectral bins
set to 4 or more. The optimal number of bins depends on numerous factors, such as
- If Mitsuba is to be run in parallel mode across a number of network compute resources
- The complexity of the scenes to be rendered
- The number of spectral bins that are actually required for the problem at hand

Mitsuba has a number of "integrators", being the plugins that actually implement different rendering schemes.
The path tracer (``path``) is the integrator to be selected for general purposes, where there is direct and indirect
illumination of the scene. For faster renders, the direct illumination integrator (``direct``) provides good quality
renders without indirect illumination components. For the path tracer, the Hammersley QMC sampler is preferred, with
as many as 256 samples per pixel or more to reduce monte carlo noise.