Installation and Requirements
=============================
 ``MORTICIA`` has been developed largely in Python 2.7 and has not yet been tested in Python 3.X.
 A working installation of `libRadtran <http://www.libradtran.org>`_ is required to compute radiant environment
 maps and atmospheric transmittance. In the Monte carlo/statistical mode of operation, a compute cluster
 is generally required to achieve adequate sampling in a reasonable time. Parallel computation is performed
 using the `ipyparallel <https://ipyparallel.readthedocs.org/en/latest/>`_ package, which works in
 IPython/Jupyter notebooks as well as other Python launch modes. In principle, with ``ipyparallel`` it is possible
 to run the Python worker engines on any network-accessible compute resource. This means that it is possible to
 utilise a Linux-based compute cluster with a libRadtran installation in the background, while all foreground
 work is performed on a Windows machine, for example. However, this setup may require
 significant setup to implement. All the required information to do this is provided in the ``ipyparallel``
 documentation. Examples are provided in the `MORTICIA` notebooks in the `nbMORITICIA repository on GitHub <https://github.com/derekjgriffith/nbMORTICIA>`_.

MORTICIA Dependencies
=====================
``MORTICIA`` has been developed using the `Anaconda <https://www.continuum.io/downloads>`_ distribution from
`Continuum Analytics <https://www.continuum.io/>`_ and this is the recommended distribution for ``MORTICIA`` users.
In principle, any Python 2.7 installation that can meet the dependencies should also work.

``MORTICIA`` makes use of a variety of other Python packages. The most important of these are:

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
- ``ipyparallel`` is used for parallel computation in IPython/Jupyter notebooks as well as other Python launch modes.
- ``dill`` will be required if running libRadtran on remote machines using ``ipyparallel``
- ``paramiko`` will be required for authentication when using ``ipyparallel`` across machine boundaries

