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

Test Heading
============

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
============

MORTICIA is Monte carlo Optical Rendering for Theatre Investigations of Capability under the Influence of the
Atmosphere. It is intended to model optical spectrum imaging sensors observing targets through the atmosphere and
to provide a statistical picture of the effectiveness of such sensors deployed in particular environments or
geographical theatres.

MORTICIA can be used to:

- Model a digital image recorded by a camera (imaging optical sensor) of a particular target under specific atmospheric conditions at a given time of day.
- Model imagery from ground-based, airborne or space optical sensors.
- Estimate the range of Detection, Recognition and Identification (DRI) for imaging optical sensors operating in the wavelength bands from the ultraviolet (~ 300 nm wavelength) to the thermal infrared (~ 14 microns wavelength).


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







