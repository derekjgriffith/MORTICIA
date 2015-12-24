__author__ = 'DGriffith'
"""
.. module:: optics
    :platform: Windows, Unix
    :synopsis: The optics module includes all code related to imaging optics as spatial and spectral filters. It also
               includes everything related to light propagation within such imaging optics. It does not include the
               atmospheric radiative transfer code. Functions related to the optical characteristics of the human
               eye are included in this module.
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




