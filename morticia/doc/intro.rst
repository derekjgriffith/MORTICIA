Introduction
============

``MORTICIA`` is Monte carlo Optical Rendering for Theatre Investigations of Capability under the Influence of the
Atmosphere. It is intended to model optical spectrum imaging sensors observing targets through the atmosphere and
to provide a statistical picture of the effectiveness of such sensors deployed in particular environments or
geographical theatres.

``MORTICIA`` can be used to:

- Model a digital image recorded by a camera (imaging optical sensor) of a particular target under specific atmospheric conditions at a given time of day.
- Model imagery from ground-based, airborne or space optical sensors.
- Estimate the range of Detection, Recognition and Identification (DRI) for imaging optical sensors operating in the wavelength bands from the
  ultraviolet (~ 300 nm wavelength) to the thermal infrared (~ 14 microns wavelength).

The most important top-level classes in ``MORTICIA`` are Sensor, Scene and Snap.

A Sensor comprises an imaging camera on a platform. A Scene comprises at least 2 scene elements (a target and a
background). A Snap comprises a single integration period of the camera (snapshot), viewing the scene in a specific geometry
(including range) under a specific set of atmospheric and illumination conditions. Illumination can be natural
(sun, moon, stars) or artificial such as using a laser illuminator.