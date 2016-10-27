# MORTICIA
Note : ``MORTICIA`` is in the very early stages of development. The following description is based on the desired result.
## Introduction
``MORTICIA`` is Monte-carlo Optical Rendering for Theatre Investigations of Capability under the Influence of the Atmosphere, 
an open-source optical remote sensing and surveillance modelling framework written mainly in Python.

  The theatre of operations is a geographical region defined by a polygon in which an optical remote sensing or
surveillance system are deployed or to be deployed. Besides the engineering characteristics of the optical remote sensing or 
surveillance system, the atmosphere plays a key role in determining system effectiveness. In the case of a 
surveillance system, effectiveness may be quantified in terms of the range at which targets can be 
Detected, Recognised or Identified (DRI). 
  
  Some of the key capabilities of ``MORTICIA`` are:

 + Physically accurate rendering of targets of interest in atmospheric conditions that are statistically representative 
  of the theatre of operations.
 + Statistically valid Monte-carlo simulations of the Detect, Recognition and Identification (DRI) ranges of such targets
 in the theatre of operations.
 
## Theatre Climatology
 The climatology within the theatre must be taken into account in order to obtain statistically valid results. Some of the
 most important elements are:
 
 + Cloud climatology
 + Aerosol climatology
 + Surface reflectance climatology (e.g. vegetation state)
 + Atmospheric turbulence climatology
 
The atmospheric turbulence climatology is especially important for high resolution, long range optical 
surveillance systems deployed within the Atmospheric Boundary Layer (ABL). High resolution space sensors are
impacted less by turbulence in general.
 
 In addition, there are certain fixed geometrical aspects relevant to surveillance that are defined by selection of the
 theatre. These include solar/lunar azimuth and elevation statistics and surface topography.
 
## Outputs
 Typical output of ``MORTICIA`` takes the form of a database of results that define target appearance and DRI, having chosen
 a sensor position and a target position within the theatre at a particular time of day and day of year. The atmospheric 
 conditions can be chosen explicitly or generated randomly on the basis of the theatre climatology.
 
 The Python package ``pandas`` is the preferred environment for mining the results.
 
## Tools
 ``MORTICIA`` tools include utilities for building 
 + theatre boundaries and climatologies,
 + theatre topography
 + target geometry and optical characteristics
 + sensor characteristics
 + radiant environment maps (using the libRadtran radiative transfer code)
 
## Documentation
 Code documentation for ``MORTICIA`` is generated using [Sphinx](http://sphinx-doc.org).
 
## Installation and Requirements
 ``MORTICIA`` has been developed largely in Python 2.7 and has not yet been tested in Python 3.X.
 A working installation of [libRadtran](http://www.libradtran.org) is required to compute radiant environment
 maps and atmospheric transmittance. In the Monte carlo/statistical mode of operation, a compute cluster
 is generally required to achieve adequate sampling in a reasonable time. Parallel computation is performed
 using the [`ipyparallel`](https://ipyparallel.readthedocs.org/en/latest/) package, which works in 
 IPython/Jupyter notebooks as well as other Python launch modes.
 
``MORTICIA`` has been developed using the [Anaconda](https://www.continuum.io/downloads) distribution from
[Continuum Analytics](https://www.continuum.io/) and this is the recommended distribution for ``MORTICIA`` users.
In principle, any Python 2.7 installation that can meet the dependencies should also work. 
 
 
## Repository
 The master repository for ``MORTICIA`` is publicly hosted on [GitHub](http://www.github.org) at 
 [https://github.com/derekjgriffith/MORTICIA](https://github.com/derekjgriffith/MORTICIA).
 
 
 
 

