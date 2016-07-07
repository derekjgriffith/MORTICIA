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
 documentation. Examples are provided in the ``MORTICIA`` notebooks in the `nbMORITICIA repository on GitHub <https://github.com/derekjgriffith/nbMORTICIA>`_.

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
- ``astropy`` and ``pyephem`` : for various astronomy-related calculations (sun, moon position). ``astropy`` may be part
  of the Anaconda distribution, which is continually growing.
- ``easygui`` is used for simple file/open and save dialogs.
- ``ipyparallel`` is used for parallel computation in IPython/Jupyter notebooks as well as other Python launch modes.
- ``dill`` will be required if running libRadtran on remote machines using ``ipyparallel``
- ``paramiko`` will be required for authentication when using ``ipyparallel`` across machine boundaries
- ``openexr`` for writing environment maps to OpenEXR files

Obtatining working installations on Windows is sometimes problematic. Unofficial Python Wheel dsitributions for
Windows can be obtained from `Christoph Gohlke <http://www.kaij.org/blog/?p=123>`_.

For example, the ``OpenEXR`` package may have to be downloaded as a Python Wheel and installed from the command-line
with ``pip``::

    C:\Ananconda> pip install OpenEXR-1.2.0-cp27-none-win_amd64.whl



Setup of ``ipyparallel``
========================
In order to run Python processes for ``MORTICIA`` in parallel on a compute cluster, the ``ipyparallel`` package is required. The
following general steps must be followed:

- Install Git, Anaconda and other ``MORTICIA`` dependencies including ``ipyparallel`` on any machines (nodes) that
  will be part of the compute cluster resource. Also install `libRadtran <http://www.libradtran.org>`_  if required on the cluster.
- Clone the ``MORTICIA`` repository to the compute nodes and point Python to the relevant installation using a .pth file. For
  example, if the ``MORTICIA`` repo has been cloned to `~/GitHub/MORTICIA` and Anaconda is installed in your home directory
  in `~/anaconda` then create a .pth file called
  `~/anaconda/lib/python2.7/site-packages/morticia.pth` which provides the full path to the ``MORTICIA`` repo clone.
- Create a working directory on each machine in the cluster. The compute engines will be started in this working
  directory. If using ``libRadtran``, the work directoy should be a sub-directory of the directory containing the
  `libRadtran` `data` directory unless other arrangements have been made for `libRadtran` to locate relevant input
   data.
- Run the `ipcontroller --ip=*` command in the working directory to start the cluster controller process listening
  on all interfaces. This should only be done if your cluster resides on a trusted network. Consult the documentation
  at `ipython.org <https://ipython.org/ipython-doc/2/parallel/parallel_process.html>`_ for further details.
  The `ipcontroller` process will announce writing of .json files which are needed for clients and compute engines to
  connect to the controller.
  Make a note of the full pathnames of these .json files. They are typically in `~/.ipython/profile_default/security` with
  names `ipcontroller-client.json` and `ipcontroller-engine.json`.
- Copy the `ipcontroller-engine.json` file to a similar location on the machine that will host the compute engines. If
  this is the same machine on which the controller itself is running, then this step may not be required.
- Copy the `ipcontroller-client.json` file to the client that will require compute resources in the
  ``IPYTHONDIR/profile_default/security` directory. On Windows this is typically `C:\Users\myself\.ipython\profile_default\`.
  Additional profile directories can be created as required.

