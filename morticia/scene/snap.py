__author__ = 'DGriffith'

"""

.. module:: snap
    :platform: Linux/Windows
    :synopsis: Module providing for imaging snapshots of targets using optical imaging sensors.

This module provides support for taking a snapshot of a scene. A snapshot comprises a single exposure
of a particular sensor, pointing into a specific scene under particular viewing geometry and
atmospheric conditions including location and time of day/year.

The CZML language is used to describe scenarios. The library at https://github.com/cleder/czml
is used to read and write CZML files in Python.

Scenario presentation is performed using the Cesium project https://cesium.com/index.html.
"""

class Snap(object):
    """
    Class encapsulating a "snapshot" of a scene. It includes information on the scenario in which the
    snapshot is to be taken, including the sensor information, target(s) and atmospheric conditions.

    There are multiple ways in which snapshots can be instantiated.

    """
    def __init__(self, czml_doc, czml_filename=None, sensor_id, ):
        """
        Construct a snapshot
        :param czml_filename:
        :return:
        """