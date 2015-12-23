# Create a global units registry
from pint import UnitRegistry
ureg = UnitRegistry()
Q_= ureg.Quantity
def U_(units):
    """ Returns a pint.Quantity of magnitude one with the given units

    :param units: A string providing a unit of measure in the pint unit registry
    :return: A pint.Quantity object with a magnitude of 1 and the specified units.
    """
    return Q_(1.0, units)