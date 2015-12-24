Sensor Package
===============
The sensor package provides modelling of complete sensors, typically comprising an optical objective lens and
a focal plane detector as a minimum (a camera). A sensor may also include a display and a human observer. The display
and human observer are required in order to calculate ranges and probabilities of target detection, recognitions and
identification (DRI).

Capability is also provided to model direct view systems comprising a telescope or binocular and human observer, without
any electronics in the sensor chain.

Optics Module
-------------
The optics module provides modelling of basic objective lenses used by imaging systems. For the most part, these are
optical systems that are near to diffraction-limited, with centred circular entrance pupils, possibly with a centred
circular obscuration (as in a Cassegrain or Newtonian telescope).

.. automodule:: optics
   :members:

Electronics Module
------------------

The electronics module provides modelling of all electronic componenents in the imaging chain. Inclusion
of any electronic component into a sensor results in an "electro-optical" sensor. Electronic components include

- Solid state electro-optical detectors, including focal plane array (FPA) detectors
- Image intensifier tubes
- Image processing components
- Display systems

