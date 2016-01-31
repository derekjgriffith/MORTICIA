Packages
========
This is a sub-menu for all the packages

.. _sensor-package:

Sensor Package
--------------
The sensor package provides modelling of complete sensors, typically comprising an optical objective lens and
a focal plane detector as a minimum (a camera). A sensor may also include a display and a human observer. The display
and human observer are required in order to calculate ranges and probabilities of target detection, recognitions and
identification (DRI).

Capability is also provided to model direct view systems comprising a telescope or binocular and human observer, without
any electronics in the sensor chain.

.. _optics-module:

Optics Module
~~~~~~~~~~~~~
The optics module provides modelling of basic objective lenses used by imaging systems. For the most part, these are
optical systems that are near to diffraction-limited, with centred circular entrance pupils, possibly with a centred
circular obscuration (as in a Cassegrain or Newtonian telescope).

.. automodule:: optics
   :members:

.. _electronics-module:

Electronics Module
~~~~~~~~~~~~~~~~~~

The electronics module provides modelling of all electronic components in the imaging chain. Inclusion
of any electronic component into a sensor results in an "electro-optical" sensor. Electronic components include

- Solid state electro-optical detectors, including focal plane array (FPA) detectors
- Image intensifier tubes
- Image processing components
- Display systems

.. automodule:: electro
   :members:
   
.. _radiometry-package:

Radiometry Package
------------------
The radiometry package provides all MORTICIA functionality relating to radiometry and atmospheric radiative transfer.
Atmospheric radiative transfer is performed by the `libRadtran <http://www.libradtran.org>`_ suite of tools.


.. _librad-module:

LibRadTran Module
~~~~~~~~~~~~~~~~~~

The LibRadTran module concerns itself with calculations for the radiant environment map. It solves the radiative transfer equation for a complex system.

- Item 1
- Item 2
- Item 3
- Item 4

.. automodule:: librad
   :members:
   

   
.. _scene-package:
   
Scene Package
------------------
The Scene package provides all MORTICIA functionality relating to radiometry and atmospheric radiative transfer.
Atmospheric radiative transfer is performed by the `libRadtran <http://www.libradtran.org>`_ suite of tools.


.. _target-module:

Target Module
~~~~~~~~~~~~~~~~~~

The Target module concerns itself with calculations for rendering a target

- Greyscale Beach Target
- Target 2
- Target 3

.. automodule:: target
   :members: