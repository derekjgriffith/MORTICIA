# Create a global units registry
# All morticia packages/modules that use the unit registry should import the global copy using
# from .. import ureg, Q_, U_
from pint import UnitRegistry
ureg = UnitRegistry()
Q_= ureg.Quantity
def U_(units):
    """ Returns a pint.Quantity of magnitude one with the given units

    :param units: A string providing a unit of measure in the pint unit registry
    :return: A pint.Quantity object with a magnitude of 1 and the specified units.
    """
    return Q_(1.0, units)