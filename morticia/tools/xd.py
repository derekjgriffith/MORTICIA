__author__ = 'DGriffith'

# Functions related to interpolation of xray.DataArray and other utilities
# The following function takes two DataArray objects and interpolates
# both onto a merged set of coordinate points.

# Some sort of global data dictionary (CF compliant ?), including short names and synonyms
# Could import some or all CF definitions from XML file.

import numpy as np
import xray
from .. import ureg, Q_, U_

def xd_harmonise_interp(dar_tup):
    """ Perform linear interpolation on merged set of axis points for two or more xray DataArray objects.
    This function can be used to prepare (harmonise) multiple xray.DataArray objects for multiplication or addition
    on a common set of coordinate axis points by linearly interpolating all DataArray objects onto the same
    set of points, obtained by merging and sorting the points from all input DataArray objects.

    The DataArry objects provided.    The scipy linear grid interpolator is used for this purpose. See:
    scipy.interpolate.RegularGridInterpolator
    :param dar_tup:
    :return: Tuple of xray.DataArray objects with merged and linearly interpolated values in all axes.
    Only unique values in the interpolation axis are used.

    """
    # TODO enforce compatible attributes or not ? What attributes in returned object

    # Accumulate the index values from each of the given arrays, for each of the axes in the first array
    index_vals = {}  # dictionary of index coordinates for each axis

    for dar in dar_tup:
        for axis in dar.dims:
            pass  # accumulate dictionary for all dimensions in the entire collection of DataArrays

    # get the unique values in increasing numerical order using np.unique for each axis found in the whole set

    # interpolate each of the DataArray objects onto the new grid (for whatever axes it does have)

    # There may be axes not present in a specific DataArray. These are omitted for that DataArray and
    # simply allowed to broadcast hen the


def check_convert_units(value_with_units, preferred_units):
    """ Check the units of a quantity and convert to preferred units using Python `pint`

    :param value_with_units: A list with a numeric value or numpy array in the first position and a string
        providing units in the second position. The unit string must be recognisable by the Python `pint` package.
    :param preferred_units: A string expressing the units to which `pint` should convert the scalar
    :return: Value expressed in the preferred units
    """

    # Use pint to convert
    value = Q_(np.asarray(value_with_units[0], dtype=np.float64), value_with_units[1])  # Will blow up if units not recognised
    value = value.to(preferred_units)
    return value.magnitude


def xd_check_convert_units(xd, axis_name, preferred_units):
    """ Check and convert units for one or more axes of an `xray.DataArray`

    :param xd: An xray.DataArray object having an axis called `axis_name` and a value in the `attrs` dictionary
    :param preferred_units: A string providing the preferred units that can be passed to `pint`
    :return: A xray.DataArray, inwhich the values in the named axis have been converted to the preferred units
        The `axis_name_units` field is also updated.
    """

    # Create a pint.Quantity object using the data from the named array
    Q_values = Q_(xd[axis_name].data, xd.attrs[axis_name + '_units'])
    Q_values = Q_values.to(preferred_units)
    xd[axis_name] = Q_values.magnitude
    xd.attrs[axis_name + '_units'] = preferred_units



