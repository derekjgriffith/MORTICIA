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
- ``easygui`` is used for simple file/open and save dialogs.
- ``xray`` is used for multi-dimensional data hypercubes with named axes.