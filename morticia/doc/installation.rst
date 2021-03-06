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
 documentation. Examples are provided in the ``MORTICIA`` notebooks in the
 `nbMORITICIA repository on GitHub <https://github.com/derekjgriffith/nbMORTICIA>`_.

 The recommended Linux platform is Debian 9 (stretch) or Ubuntu based on Debian stetch (16.04 or later).

MORTICIA Dependencies
=====================
``MORTICIA`` has been developed using the `Anaconda <https://www.continuum.io/downloads>`_ distribution from
`Continuum Analytics <https://www.continuum.io/>`_ and this is the recommended distribution for ``MORTICIA`` users.
In principle, any Python 2.7 installation that can meet the dependencies should also work.

``MORTICIA`` makes use of a variety of other Python packages. Many `MORTICIA` dependencies are included in the
Anaconda package. The most important of are:

- ``numpy`` and ``scipy`` : The mainstay of technical computing in Python.
- ``matplotlib`` : Plotting.
- ``xarray`` : A package for handling N-dimensional data arrays and datasets with named axes. Also for reading and
writing data from netCDF format files.
- ``pandas`` : For handling time-series and data in a tabular/relational view.
- ``lxml`` : For reading, writing, querying and transforming files in the `XML` format (eXtensible Markup Language).
- ``pint`` : for handling unit checking and conversion. Note that using ``pint`` for all numerical operations is not
  recommended and not done in ``MORTICIA``. Rather, a ``pint`` Quantity object can be carried as metadata to an
  xarray.DataArray object and checked (for example using Python ``assert``) at object construction time or judiciously
  so that there is no significant performance impact of any unit checking.
- ``astropy`` and ``pyephem`` : for various astronomy-related calculations (sun, moon position). ``astropy`` may be part
  of the Anaconda distribution, which is continually growing.
- ``easygui`` is used for simple file/open and save dialogs.
- ``ipyparallel`` is used for parallel computation in IPython/Jupyter notebooks as well as other Python launch modes.
- ``dill`` will be required if running libRadtran on remote machines using ``ipyparallel``
- ``paramiko`` will be required for authentication when using ``ipyparallel`` across machine boundaries
- ``openexr`` for writing environment maps to OpenEXR files

For physically based and detailed optical rendering, ``MORTICIA`` interfaces to the
`Mitsuba <http://www.mitsuba-renderer.org/>_`  rendering system.


Other useful auxiliary applications include:

- `mrViewer <http://mrviewer.sourceforge.net/>`_ for viewing floating point or multi-channel OpenEXR files.
- `Blender <https://www.blender.org/>`_ for 3D model import, editing, texturing and export to Collada format for
import into ``Mitsuba``.
- `MeshLab <http://www.meshlab.net/>`_ for import, visualisation, analysis, repair and format conversion of geometry
meshes
for import into `Blender` or directly into `Mitsuba`.

Obtaining working installations for all of these packages on Windows is sometimes problematic. Unofficial Python Wheel
distributions for
Windows can be obtained from `Christoph Gohlke <http://www.kaij.org/blog/?p=123>`_.

For example, the ``OpenEXR`` package may have to be downloaded as a Python Wheel and installed from the command-line
with ``pip``::

    C:\Ananconda> pip install OpenEXR-1.2.0-cp27-none-win_amd64.whl

Environment maps written in the OpenEXR format are used for integration with the Mitsuba rendering system.

Git
----

The code management and revision control system `Git` is used for source control in the `MORTICIA` and several
dependent projects e.g. Mitsuba.

Install Git on the relevant platform(s) using `sudo apt-get install git` on Debian or Ubuntu or by downloading an
installer from the `Git download page <https://git-scm.com/downloads>`_.

Mitsuba
-------

For physically-based rendering of target scenes, the supported rendering system is `Mitsuba <http://www
.mitsuba-renderer.org/>`_.

Mitsuba can be obtained in several forms, including simple installers for Windows at the `Mitsuba download
page <https://www.mitsuba-renderer.org/download.html/>`_. For serious multi-spectral or hyperspectral work, Mitsuba
must be compiled from source. The source code is also available from the Mitsuba download page, but the `Mitsuba Github
repository <https://github.com/mitsuba-renderer/mitsuba>`_ may have more recent code and bug fixes. Subscribe to
Watch on the Mitsuba Github repository to stay informed about the latest activity on Mistuba development. If
 the Git repository is used, it is important to obtain also the latest Mitsuba manual, which may have to be built
 from the `TeX` source. The latest manual also contains the latest build instruction.

For multispectral or hyperspectral work, Mitsuba must be compiled from source in the spectral mode with the number of
spectral bins set to 4 or more. The Mitsuba manual contains detailed instructions for compilation of Mitsuba and how
to set the number of spectral bins (SPECTRAL_SAMPLES compilation flag). The optimal number of bins depends on numerous
factors,
such as
- If Mitsuba is to be run in parallel mode across a number of network compute resources
- The complexity of the scenes to be rendered
- The number of spectral bins that are actually required for the problem at hand

When compiling Mitsuba on a remote compute server, it is best to add the MTS_GUI_SOFTWARE_FALLBACK compilation
 flag. See the Mitsuba manual for further details.

Mitsuba has a number of "integrators", being the plugins that actually implement different rendering schemes.
The path tracer (``path``) is the integrator to be selected for general purposes, where there is direct and indirect
illumination of the scene. For faster renders, the direct illumination integrator (``direct``) provides good quality
renders without indirect illumination components. For the path tracer, the Hammersley QMC sampler is preferred, with
as many as 256 samples per pixel or more to reduce monte carlo noise.

The Collada Document Object Model (DOM) allows for Mitsuba to make use of 3D model geometry in the Collada `.dae`
file format. Mitsuba is therefore preferably compiled including the Collada DOM as per the instructions in the Mitsuba
 manual. This is not mandatory provided that 3D object models are available in file formats that are natively
 supported in Mitsuba (`.serial` format).

See the Mitsuba documentation for further details.

A specific limitation with Mitsuba is that compilations with different numbers of spectral bins are not compatible
with one another. It is useful to have an RGB version of Mitsuba available to compute distance maps of a scene. This
is best done on the local host, which should have an RGB version of Mitsuba installed (and will therefore not bw able
 to cooperate using `mtssrv` with remote compute nodes).

Issues with Mitsuba on Debian 9
===============================

Currently (Jan 2018), Debian 9 (stretch) is not yet officially supported by Mituba.

Debian 9 (stretch) is the current preferred platform for remote compute nodes whe using `MORTICA`. Note the following
with respect to Debian 9:
- Use the Github distribution of Mitsuba and compile from source
- Do not install Qt or prevent Mitsuba from attempting to compile the GUI (`mtsgui`) by manually setting the hasQt
variable to False in the file build/SConscript.configure. Attempting to compile the GUI results in Mitusba build
errors on Debian stretch.
- Do not install the Collada DOM. Currently it causes Mitsuba build errors on Debian stretch. Building the Collada DOM
from source may work. However, it is best to convert Collada models to the much more efficient Mitusba `.serial` format
before using them in heavy modelling work. This can be done using a local install of RGB Mistuba and `mtsgui`.
- Build and install Mitsuba before Anaconda.
- Build OpenEXR from source and perform a system-wide install (IlmBase and OpenEXR). The stretch repository version
does not appear to work with Mitsuba.
- Ubuntu versions from 16.04 are based on Debian stretch and may therefore be suitable for `MORTICIA` remote compute
node purposes.

Cesium and CZML
===============
`Cesium <http://cesiumjs.org>`_ is a javascript/browser-based system for 4D (XYZT) simulation and display. An example
 of the Cesium viewer can be seen at the `CesiumJS website <https://cesiumjs.org/Cesium/Build/Apps/CesiumViewer/index.html>`_



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
   data. The recommended name of this directory is `work`.
- Open a detachable terminal using the `screen` command or using the `byobu` package. Use one detachable screen in
  which to start the `ipcontroller` and another in which to run the engines. Rename these screen sessions to make them
  easier to locate later on, using `ctrl-a :sessionname ipcontroller`.
- Run the `ipcontroller --ip=*` command in the working directory to start the cluster controller process listening
  on all interfaces. This should only be done if your cluster resides on a trusted network. Consult the documentation
  at `ipython.org <https://ipython.org/ipython-doc/2/parallel/parallel_process.html>`_ for further details.
  The `ipcontroller` process will announce writing of .json files which are needed for clients and compute engines to
  connect to the controller.
  Make a note of the full pathnames of these .json files. They are typically in `~/.ipython/profile_default/security` with
  names `ipcontroller-client.json` and `ipcontroller-engine.json`.
- Detach from the `ipcontroller` screen session.
- Copy the `ipcontroller-engine.json` file to a similar location on the machine that will host the compute engines. If
  this is the same machine on which the controller itself is running, then this step may not be required.
- Copy the `ipcontroller-client.json` file to the client that will require compute resources in the
  ``IPYTHONDIR/profile_default/security` directory. On Windows this is typically `C:\Users\(your_username)\.ipython\profile_default\`.
  Additional profile directories can be created as required.
- If the controller is restarted for any reason, the above `.json` files will be overwritten and the above file copy
  operations must be repeated. This can be avoided by using `ipcontroller --reuse`.
- The compute engines can now be started on the relevant machines in the desired working directory using
  the `ipengine` command. Preferably start another detachable `screen` or `byobu` session for this and rename the
  session to `ipengine`.
- Note that it is generally very important to ensure that the same version of all dependent Python packages is
  running on all compute nodes. Once all dependencies have been installed, make sure that all packages are updated
  or upgraded, in particular the ``xarray`` package (`pip install --upgrade xarray`).
- On Windows, it may be necessary to run the command window as Administrator to get the necessary privileges for
  package installation and upgrading.

It is very important to keep the ``MORTICIA`` code the same on all platforms in use. Git pull the code and restart the
compute engines on the compute nodes if the ``MORTICIA`` codebase is altered. Also restart on the host. A common
symptom of code that is out of sync when using `ipyparallel` is a PicklingError exception.

Setup of Dask Distributed
=========================
The ``dask`` and ``distributed`` packages provide distributed parallel processing similar to that provided by
``ipyparallel``, but without security of any kind. Because ``dask.distributed`` is so easy to use, it is the
preferred option for trusted networks. See the example Jupyter notebooks on calculation of radiant environment maps
using libRadtran.



MORTICIA Development
====================
Contributing to MORTICIA or the MORTICIA notebooks requires setup of a development environment using ``conda``.
Depending on whether development is done in the Python 2.7 or a Python 3.3 context, this may entail creating a
development environment called mordevpy27 or mordevpy33. If the development environment is to use the same Python
as the Anaconda root environment, this can be done as a simple clone. Once a full installation of Anaconda has
been completed, an Anaconda or normal terminal can be opened and the development environment cloned from the root
using::

    > conda create --name mordevpy27 --clone root

Once the clone has completed, the environment can be activated with::

    > activate mordevpy27

in Windows or::

    > source activate mordevpy27

in Linux.

Not all packages required by `MORTICIA` or `nbMORTICIA` are included with Anaconda. These will have to be installed
manually using `conda` or `pip`. Missing packages may include `paramiko`, `pint`, '`easygui`, `dill`, `ipyparallel` and
 `xarray`. If the development environment is not cloned from root, it will be necessary to install many more
 packages, including basics
such as `numpy`.

Keeping Installations Synchronised
==================================
`MORTICIA` generally requires multiple platforms, including a remote compute server/cluster. It is necessary to keep
all components up-to-date on all platforms.
