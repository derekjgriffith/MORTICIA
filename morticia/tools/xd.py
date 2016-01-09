__author__ = 'DGriffith'

# Functions related to interpolation of xray.DataArray and other utilities
# The following function takes two DataArray objects and interpolates
# both onto a merged set of coordinate points.

# Some sort of global data dictionary (CF compliant ?), including short names and synonyms
# Could import some or all CF definitions from XML file.

import numpy as np
import xray
from .. import ureg, Q_, U_
from scipy.interpolate import RegularGridInterpolator, interpn
import warnings
from operator import mul, add
from ..moglo import *  # Import morticia global vocab, exceptions etc.

def xd_identity(np_vector, axis_name, attrs=None):
    """ Create an identity xray.DataArray. That is, a DataArray vector in which both the values and axis
        coordinates are identical.

    :param np_vector: Vector of numeric data
    :param axis_name: Name for the axis - must be in vocabulary defined in moglo.py
    :param attrs: Dictionary of additional attributes to attach to the DataArray
    :return:
    """
    if axis_name in long_name:
        the_long_name = long_name[axis_name]
    else:
        warnings.warn('Unknown axis name ' + axis_name + ' encountered in xd_identity creation.')
    if axis_name in default_units:
        the_units = default_units[axis_name]
    else:
        the_units = ''
    if attrs is not None:
        the_attrs = attrs
    else:
        the_attrs = {}
    the_attrs.update({'long_name': the_long_name})
    the_attrs.update({'units': the_units})
    return xray.DataArray(np_vector, [(axis_name, np_vector)], name=axis_name, attrs=the_attrs)


def xd_harmonise_interp(xd_list):
    """ Perform linear interpolation on merged set of axis points for two or more xray DataArray objects.
        This function can be used to prepare (harmonise) multiple xray.DataArray objects for multiplication or addition
        on a common set of coordinate axis points by linearly interpolating all DataArray objects onto the same
        set of points, obtained by merging and sorting the points from all input DataArray objects.

    The DataArry objects provided.    The scipy linear grid interpolator is used for this purpose. See:
    scipy.interpolate.RegularGridInterpolator
    :param xd_list:
    :return: Tuple of xray.DataArray objects with merged and linearly interpolated values in all axes.
    Only unique values in the interpolation axis are used.

    """
    # TODO : enforce compatible attributes or not ? What attributes in returned object ?
    # TODO : Ignore axes that have non-numeric coordinates e.g. xdarray['axisname'].dtype.char in 'SUa',
    # TODO : which detects dtypes that are string, or xdarray['axisname'].dtype.kind in 'fc' (float or complex)
    # TODO : alternatively require that axis coordinates are always numeric, say with a list of labels as attrs
    # TODO : What about interpolation on times axes
    # TODO : Need to expand on extrapolation (and possibly also single-axis interpolation) schemes
    # Accumulate the index values from each of the given arrays, for each of the axes in the first array
    index_vals = {}  # dictionary of index coordinates for each axis
    index_float = {}  # determine if the index kind is a floating point type (complex included)
    #metadata = {}

    for xd_arr in xd_list:
        for axis in xd_arr.dims:
            # accumulate dictionary for all dimensions in the entire collection of DataArrays
            if not axis in index_vals:
                index_vals[axis] = xd_arr[axis]
            else:
                index_vals[axis] = np.hstack((index_vals[axis], xd_arr[axis]))
            # also accumulate the attributes (metadata)
            # metadata.update(xd_arr.attrs)
    # get the unique values in increasing numerical order using np.unique for each axis found in the whole set
    for axis in index_vals:
        index_vals[axis] = np.unique(index_vals[axis])
        index_float[axis] = index_vals[axis].dtype.kind in 'fc'
    # interpolate each of the DataArray objects onto the new grid (for whatever axes it does have)
    xd_return_list = []
    for xd_arr in xd_list:
        # Create the linear interpolator
        interpolator = RegularGridInterpolator([xd_arr[axis].values for axis in xd_arr.dims],
                                               xd_arr.values,
                                               method='linear', bounds_error=False, fill_value=0.0)
        merged_coordinates = np.meshgrid(*[index_vals[axis] for axis in xd_arr.dims],
                                               indexing='ij')
        interp_vals = interpolator(tuple(merged_coordinates))
        # reconstruct the xray.DataArray with interpolated data
        xd_arr_interp = xray.DataArray(interp_vals, [(axis, index_vals[axis]) for axis in xd_arr.dims],
                                       name=xd_arr.name, attrs=xd_arr.attrs)
        xd_arr_interp.attrs = xd_arr.attrs  # transfer the attributes verbatim
        xd_return_list.append(xd_arr_interp)

        # There may be axes not present in a specific DataArray. These are omitted for that DataArray and
        # simply allowed to broadcast when performing operations with other DataArrays
    return xd_return_list

def xd_harmonised_product(xd_list):
    """ Compute the harmonised product of a number of N-dimensional data arrays.
        The DataArrays are interpolated onto a common set of coordinates and then the product of the DataArrays
        is computer, returning a single DataArray with merged attributes. Unit mismatches are flagged with warnings.

    :param xd_list: List/tuple of xray.DataArray objects to be multiplied
    :return: Product of xray.DataArray objects with merged attributes
    :except UnitMismatch, MissingUnits:
    """
    # TODO : This function to be checked to correct "var_units" mistake

    #main_attrs = {}  # Will accumulate all main attributes here - not sure what to do with units ?
    unit_dict = {}  #  Dictionary of units
    axis_attrs = {}  # Dictionary of axis attribute dictionaries
    # Check units and merge metadata
    # have to merge attributes for main data and all axes individually
    for xd_arr in xd_list:
        #main_attrs.update(xd_arr.attrs)
        for axis in xd_arr.dims:
            if axis in axis_attrs:
                axis_attrs[axis].update(xd_arr[axis].attrs)  # Accumulate the attributes for each axis
            else:
                axis_attrs[axis] = xd_arr[axis].attrs
            if not axis in unit_dict:
                if 'units' in xd_arr[axis].attrs:
                    unit_dict[axis] = xd_arr[axis].attrs['units']
                else:
                    raise MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis)
            elif ('units' in xd_arr[axis].attrs['units']) and (unit_dict[axis] != xd_arr[axis].attrs['units']):
                # TODO : Consider throwing a unit mismatch error, or converting to desired units with pint
                warnings.warn('Unit mismatch found when taking xray.DataArray harmonised product.')
                raise UnitMismatch('Unit mismatch encountered for ' + xd_arr.name + ' on axis ' + axis)
            else:
                raise MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis)
    xd_factors = xd_harmonise_interp(xd_list)
    xd_product = reduce(mul, xd_factors)  # take the product by reducing the list using the mul operator
    #xd_product.attrs = main_attrs
    for axis in xd_product:  # Put the merged attributes into each of the axes
        xd_product[axis].attrs = axis_attrs[axis]
    return xd_product


def check_convert_units(value_with_units, preferred_units):
    """ Check the units of a quantity and convert to preferred units using Python `pint`

    :param value_with_units: A list with a numeric value or numpy array in the first position and a string
        providing units in the second position. The unit string must be recognisable by the Python `pint` package.
    :param preferred_units: A string expressing the units to which `pint` should convert the scalar
    :return: Value expressed in the preferred units
    """

    # Use pint to convert
    try:
        value = Q_(np.asarray(value_with_units[0], dtype=np.float64), value_with_units[1])  # Will blow up if units not recognised
    except TypeError:
        warnings.warn('A scalar value was supplied without units. Example of correct scalar input is [40.0, "degC"]')
    value = value.to(preferred_units)
    return value.magnitude


def xd_check_convert_units(xd, axis_name, preferred_units):
    """ Check and convert units for one or more axes of an `xray.DataArray`

    :param xd: An xray.DataArray object having an axis called `axis_name` and a value in the `attrs` dictionary
    :param axis_name: Name of the axis or data to check/convert
    :param preferred_units: A string providing the preferred units that can be passed to `pint`
    :return: A xray.DataArray, in which the values in the named axis have been converted to the preferred units
        The `axis_name_units` field is also updated.
    """

    # Create a pint.Quantity object using the data from the named array and use that to convert to
    # preferred units
    if axis_name == xd.name:  # The primary data must be converted
        Q_values = Q_(xd.data, xd.units)  # Can fetch units this way, but not set them
        Q_values = Q_values.to(preferred_units)
        xd.data = Q_values.magnitude
        xd.attrs['units'] = preferred_units
    else:  # Convert units of the named axis
        Q_values = Q_(xd[axis_name].data, xd[axis_name].units)
        Q_values = Q_values.to(preferred_units)
        xd[axis_name].data = Q_values.magnitude
        xd[axis_name].attrs['units'] = preferred_units


def xd_attrs_update(xd_list):
    """ Update long_name and units attributes of all axes in an xray.DataArray
    :param xd: Input xray.DataArray
    :return:
    """

    for xd_arr in xd_list:
        xd_arr.attrs['long_name'] = long_name[xd_arr.attrs['name']]
        for axis in xd_arr.dims:
            xd_arr[axis].attrs['long_name'] = long_name[axis]  # Will blow up if axis mnemonic name not found
            xd_arr[axis].attrs['units'] = default_units[axis]  # Likewaise if units not found for this axis name


