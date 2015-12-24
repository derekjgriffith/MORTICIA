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
These functions perform axis harmonisation and unit checking and conversion.

Units of Measure
----------------
Tracking of units of measure is not performed automatically in ``MORTICIA`` such as with the use of the pint package.
However, many functions and classes in ``MORTICIA`` expect xray.DataArray instances that provide units of measure in
the metadata. When units are provided in this way, they should be provided in manner consistent with the Pyhon units
package ``pint`` (see https://pint.readthedocs.org/ )

Commonly, when a ``MORTICIA`` function or method requires a scalar numeric input, it must be provided as a list
with magnitude and units e.g. [30, 'mm']. Unitless quantities are provided as a simple numeric magnitude.

If a variable is named ``x``, then the units for the variable can be stored as ``x_units``. This convention is
applied in the attributes of xray.DataArray objects, where, if a dimension or data is named ``var``, then attributes
are provided in a dictionary with key values ``var_units`` which contains a units string that can be provided to a
``pint Quantity`` object.

A pint global unit registry is created when `morticia` or any sub-package is imported. Any other `morticia' paackages
or modules share a single global unit registry called ureg. Convenience functions ``Q_`` for ``Quantity`` and ``U_`` for
``Quantity(1.0, unit_str)`` are also defined. Examples of usage are provided in the Jupyter notebooks/tutorials.

Logging Warnings and Exception Handling
======================
As a rule, `MORTICIA` does not use logging to files. Preferably, if any checking is performed, exceptions are thrown.
Informational messages should printed to the terminal using the `logging.info()` or `logging.debug()` calls.
Warnings that the user should take action on are provided through `warnings.warn()`. If the warning may relate to
the fact that `MORTICIA` code should be improved or debugged, it should be issued through `logging.warning()`.
`MORTICIA' internal modules do not define any logging handlers or filters. This is left to the user.
