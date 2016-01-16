MORTICIA Conventions
====================

There is a data dictionary for ``MORTICIA`` which defines commonly used variable and data names.

xray.DataArray objects are created routinely, so the constructor is abbreviated to xD.

``MORTICIA`` often deals with spectral variables, such as spectral transmission, spectral radiance etc.
More often than not, these arrays are multidimensional. The `xray` package is used in `MORTICIA` for
purposes of representing N-dimensional arrays (the `xray.DataArray` class).

However, `xray` does not (yet) support interpolation when the coordinates of the data hypercube are
not the same and two or more `xray.DataArray` objects must be added, divided or multiplied.

Utility functions that operate on xray.DataArray objects are found in the `morticia.tools.xd` package.
These functions perform axis harmonisation as well as unit checking and conversion.

The axis harmonisation functions for `xray.DataArray` objects requires that all axes have numeric coordinates.
This applies also to axes that would benefit from having text labels. Text labels would be more appropriate for
axes such as spectral channels and field orientation for MTF data. In the case of spectral channels (with axis label
`chn`) will be channel numbers. These channel numbers would normally start from 0, but in some cases, such as when
using a sub-range of correlated-k band models (e.g. `kato` or `fu` in libRadtran) the channel numbers would be a
range of integers. The convention in ``MORTICIA`` is that the text labels then be placed in an axis attribute called
`labels`.

This curently applies equally well to timeseries. For incorporation in multi-axis `xray.DataArray` objects, timestamps
must be converted to Julian date. A list of labels for such axes can and should be maintained in the `xray.DataArray`
attributes (`attrs`).

Units of Measure
----------------
Tracking of units of measure is not performed automatically in ``MORTICIA`` such as with the use of the pint package.
However, many functions and classes in ``MORTICIA`` expect xray.DataArray instances that provide units of measure in
the metadata. When units are provided in this way, they should be provided in manner consistent with the Python units
package `pint` (see https://pint.readthedocs.org/ )

Commonly, when a ``MORTICIA`` function or method requires a scalar numeric input, it must be provided as a list
with magnitude and units e.g. [30, 'mm']. Unitless quantities are provided as a simple numeric magnitude.

If a variable is named ``x``, then the units for the variable can be stored as ``units`` in the xray.DataArray
attributes dictionary (actually an OrderedDict). This units string should be recognizable to `pint` as a valid unit.
When exporting xray.DataArray objects to NetCDF files in xray.Dataset objects, this convention of using the
`units` attribute is recognised by NetCDF browsers and other NetCDF utilities. Further information on conventional
NetCDF attributes can be found at
`UCAR/Unidata <https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/metadata/DataDiscoveryAttConvention.html>`_.
See the ``MORTICIA`` notebooks for examples of creating `xray` objects with unit attributes.

The UCAR/Unidata attributes for NetCDF file elements considered as highly recommended are

- `units`, units expressed as a string - SI units preferred
- `long_name`, such as may commonly be used to label graph axes
- `standard_name`, the standard name that may be used in a specific scientific community

and for ``MORTICIA``, also

- `title`, such as may be used on the title of a plot of the data
- `labels`, such as may be used in the legends of a plot

The `standard_name` should come from the CF (Climate and Forecasting) name glossary, or the project should have
its own vocabulary of short, standard and long names for all variables. These are to be found in the
`morticia.moglo.py` module and may be expanded or updated from time to time.

A pint global unit registry is created when `morticia` or any sub-package is imported. Any other `morticia` paackages
or modules share a single global unit registry called ureg. Convenience functions ``Q_`` for ``Quantity`` and ``U_`` for
``Quantity(1.0, unit_str)`` are also defined. Examples of usage are provided in the Jupyter notebooks/tutorials.

Coordinate Systems
==================
The general observer/target topocentric, cartesian coordinate system used in `MORTICIA` is +x towards North, +y towards
the East and +z towards the zenith. Zenith angles are polar angles measured from the +z axis. This is a left-handed
system and the alternative right-handed sytem is the same, except with +z towards nadir.

The global Earth Centered Earth Fixed (ECEF) coordinate system has the north pole in the +Z direction, the prime
meridian (0 deg longitude) in the +X direction and the +Y direction at 90 deg east longitude. This is a right-handed
coordinate system.

When dealing with sightlines and camera orientations, several conventions can be used. The
`Euler angles <https://en.wikipedia.org/wiki/Euler_angles>`_ provide (relative to a given coordinate system,
typically topocentric/geodetic) the view azimuth angle (VAZ = :math:`\alpha`) and view zenith angle (VZA =
:math:`\beta`). The camera roll is then the third Euler angle :math:`\gamma`.

Logging Warnings and Exception Handling
=======================================
As a rule, `MORTICIA` does not use logging to files. Preferably, if any checking is performed, exceptions are thrown.
Informational messages should printed to the terminal using the `logging.info()` or `logging.debug()` calls.
Warnings that the user should take action on are provided through `warnings.warn()`. If the warning may relate to
the fact that `MORTICIA` code should be improved or debugged, it should be issued through `logging.warning()`.
`MORTICIA` internal modules do not define any logging handlers or filters. This is left to the user.

General Terminology and Conventions
===================================

Camera and Imager
-----------------
A `Camera` object in MORTICIA does *not* include a Lens - it is a camera body, including an FPA and a
digitisation stage (which may be fully integrated into the FPA chip - a so-called "camera on a chip").
The `Imager` class incorporates a Camera and a Lens and is capable of producing actual imagery.

Modulation Transfer Function Calculation
----------------------------------------
For image modelling purposes, it is necessary to know the full 2D MTF, whereas the normal situation is that the
MTF is described using horizontal and vertical MTF profile functions. Suppose that :math:`\eta` and :math:`\xi` are
the spatial frequencies in the horizontal and vertical directions respectively and the :math:`M\!T\!F_\eta(\eta)` and
:math:`M\!T\!F_\xi(\xi)` are the 1D MTFs in the horizontal and vertical directions. The 2D MTF is then computed as
a rotationally weighted mean of the 1D MTFs as

.. math::
    M\!T\!F(\eta,\xi)=\frac{\eta^{2}M\!T\!F_{\eta}\left(\sqrt{\eta^{2}+\xi^{2}}\right)+\xi^{2}M\!T\!F_{\xi}\left(\sqrt{\eta^{2}+\xi^{2}}\right)}{\eta^{2}+\xi^{2}},

and if the sagittal spatial frequency is defined as :math:`\rho=\sqrt{\eta^{2}+\xi^{2}}`, then

.. math::
    M\!T\!F(\eta,\xi)=\frac{\eta^{2}M\!T\!F_{\eta}\left(\rho\right)+\xi^{2}M\!T\!F_{\xi}\left(\rho\right)}{\rho^{2}}.

It is assumed here that the horizontal and vertical MTFs are symmetrical about the origin.






