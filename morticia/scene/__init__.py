__author__ = 'Ari Ramkilowan'
"""
.. module:: target
    :platform: Windows, Unix
    :synopsis: The target module includses all code used to described and render various forms of targets. Targets in this instance
                refers to an object that is under inverstigation, i.e. it is what the user is looking for, specifically. As opposed to what the user
                sees anyway
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




