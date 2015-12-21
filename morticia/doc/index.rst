.. MORTICIA documentation master file, created by
   sphinx-quickstart on Mon Apr 27 15:37:59 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MORTICIA's documentation!
====================================

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
============

``MORTICIA`` is Monte carlo Optical Rendering for Theatre Investigations of Capability under the Influence of the
Atmosphere. It is intended to model optical spectrum imaging sensors observing targets through the atmosphere and
to provide a statistical picture of the effectiveness of such sensors deployed in particular environments or
geographical theatres.

``MORTICIA`` can be used to:

- Model a digital image recorded by a camera (imaging optical sensor) of a particular target under specific atmospheric conditions at a given time of day.
- Model imagery from ground-based, airborne or space optical sensors.
- Estimate the range of Detection, Recognition and Identification (DRI) for imaging optical sensors operating in the wavelength bands from the
  ultraviolet (~ 300 nm wavelength) to the thermal infrared (~ 14 microns wavelength).

The most important top-level classes in ``MORTICIA`` are Sensor, Scene and Snap.

A Sensor comprises an imaging camera on a platform. A Scene comprises at least 2 scene elements (a target and a
background). A Snap comprises a single integration period of the camera (snapshot), viewing the scene in a specific geometry
(including range) under a specific set of atmospheric and illumination conditions. Illumination can be natural
(sun, moon, stars) or artificial such as using a laser illuminator.

MORTICIA Conventions
====================

There is a data dictionary for ``MORTICIA`` which defines commonly used variable and data names.

xray.DataArray objects are created routinely, so the constructor is abbreviated to xD.

Units of Measure
----------------
Tracking of units of measure is not performed automatically in ``MORTICIA`` such as with the use of the pint package.
However, many functions and classes in ``MORTICIA`` expect xray.DataArray instances that provide units of measure in
the metadata. When units are provided in this way, they should be provided in manner consistent with the Pyhon units
package ``pint`` (see https://pint.readthedocs.org/ )

Commonly, when a ``MORTICIA`` function or method requires a scalar numeric input, it must be provided as a list
with magnitude and units e.g. [30, 'mm']. Unitless quantities are provided as a simple numeric magnitude.

If a variable is named ``x``, then the units for the variable can be stored as ``units_x``. This convention is
applied in the attributes of xray.DataArray objects, where, if a dimension or data is named ``var``, then attributes
are provided in a dictionary with key values ``units_var`` which contains a ``pint Quantity`` object. The magnitude of
``pint Quantity`` object must be set to 1.0. Convenience functions ``Q_`` for ``Quantity`` and ``U_`` for
``Quantity(1.0, unit_str)`` are defined. Examples of usage are provided in the Jupyter notebooks/tutorials.

MORTICIA Dependencies
=====================
MORTICIA makes use of a variety of other Python packages. The most important of these are:

- ``numpy`` and ``scipy`` : The mainstay of technical computing in Python.
- ``matplotlib`` : Plotting.
- ``xray`` : A package for handling N-dimensional data arrays and datasets with named axes. Also for reading and writing
  data from netCDF format files.
- ``pandas`` : For handling time-series and data in a tabular/relational view.
- ``pint`` : for handling unit checking and conversion. Note that using ``pint`` for all numerical operations is not
  recommended and not done in ``MORTICIA``. Rather, a ``pint`` Quantity object can be carried as metadata to an
  xray.DataArray object and checked (for example using Python ``assert``) at object construction time or judiciously
  so that there is no significant performance impact of any unit checking.
- ``astropy`` and ``pyephem`` : for various astronomy-related calculations (sun, moon position).

Jupyter Notebooks
=================
Jupyter notebooks are provided to illustrate usage of the main package (``morticia``) as well as sub-packages
documented below.

Sensor Package
===============


Optics Module
-------------
The optics module provides modelling of basic objective lenses used by imaging systems. For the most part, these are
optical systems that are near to diffraction-limited, with centred circular entrance pupils, possibly with a centred
circular obscuration (as in a Cassegrain or Newtonian telescope).

.. automodule:: optics
   :members:

Radiometry Package
==================
The radiometry package provides all MORTICIA functionality relating to radiometry and atmospheric radiative transfer.
Atmospheric radiative transfer is performed by the libRadtran suite of tools


.. automodule:: librad
   :members:







