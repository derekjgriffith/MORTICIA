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

- `units`
- `long_name`
- `standard_name`

The `standard_name` should come from the CF (Climate and Forecasting) name glossary, or the project should have
its own vocabulary of short, standard and long names for all variables.

A pint global unit registry is created when `morticia` or any sub-package is imported. Any other `morticia` paackages
or modules share a single global unit registry called ureg. Convenience functions ``Q_`` for ``Quantity`` and ``U_`` for
``Quantity(1.0, unit_str)`` are also defined. Examples of usage are provided in the Jupyter notebooks/tutorials.

Logging Warnings and Exception Handling
======================
As a rule, `MORTICIA` does not use logging to files. Preferably, if any checking is performed, exceptions are thrown.
Informational messages should printed to the terminal using the `logging.info()` or `logging.debug()` calls.
Warnings that the user should take action on are provided through `warnings.warn()`. If the warning may relate to
the fact that `MORTICIA` code should be improved or debugged, it should be issued through `logging.warning()`.
`MORTICIA' internal modules do not define any logging handlers or filters. This is left to the user.
