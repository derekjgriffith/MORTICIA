__author__ = 'Ari Ramkilowan'
"""
.. module:: target
    :platform: Windows, Unix
    :synopsis: The target module includses all code used to described and render various forms of targets. Targets in this instance
    refers to an object that is under inverstigation, i.e. it is what the user is looking for, specifically. As opposed to what the user
    sees anyway
Dependencies : numpy (as np), pandas (as pd) and xray
               easygui, pint, warnings

"""

# Set default logging handler to avoid "No handler found" warnings.
''' # Currently considering using the warnings module
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())
'''
# Logging can be selectively silenced by using the method logging.Logger.setLevel() or
# disabled by setting the attribute logging.Logger.disabled to True.
from .. import moglo


import numpy as np
import pandas as pd
import xray
import warnings
import logging  # TODO set up global logging if required.

# Import units registry from parent to avoid duplicates
from morticia import ureg, Q_, U_

# Import tools related to xray DataArray and related unit checking/conversion
from morticia.tools.xd import *

class Target(object):
    """ This is a public, abstract class. All instantiations of targets will inherit and or modify
        a subset of the methods inside of the Target class
    """

    def __init__(self, name, desc, npix ):
        """Target constructor
        The target is constructed using the characteristic dimension....
        :param name: Name of target
        :param desc: Description of target
        :param npix: The number pixels at which to render target

        :return:
        """
        self.name = name
        self.desc = desc
        self.npix = npix

    def get_name(self):
        return self.name

    def get_desc(self):
        return self.desc

class GreyLambertianBeachballTarget(Target):
    """
        This class inherits the generic(abstract) Target class and adds to
        it specialised methods specific to this type of target. An instantiation
        of this class will render a greyscale beachball.
    """
    def __init__(self, diam= 1, Nsectors= 12, normal_vector= (0.0, 0.0, 0.0)):
        """

        :param name:
        :param npix:
        :return:
        """
        Target.__init__(self, 'BeachBall','Grey Spherical Lambertian',255)
        self.diam= diam
        self.Nsectors= Nsectors
        self.normal_vector= normal_vector

    def get_diam(self):
        return self.diam

    def get_Nsectors(self):
        return self.Nsectors

    def get_normal_vector(self):
        return self.normal_vector()

    def __str__(self):
        printString = """ The target is a %s %s with %i sectors, rendered with %i pixels.\n
        %.2f, %.2f, %.2f are the x, y, z componnets of the normal vector, respectively.""" \
        % (self.desc, self.name, self.Nsectors, self.npix, self.normal_vector[0], self.normal_vector[1],
           self.normal_vector[2])

        return str(printString)