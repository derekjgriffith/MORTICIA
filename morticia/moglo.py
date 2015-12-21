__author__ = 'DGriffith'

import numpy as np
import pandas as pd
import xray
import scipy.interpolate

# MORTICIA globals and utility functions

# Functions related to interpolation of xray.DataArray
# The following function takes two DataArray objects and interpolates
# both onto a merged set of points on a particular named axis

# Some sort of global data dictionary (CF compliant ?), including short names and synonyms
# Could import some or all CF definitions from XML file.

def darray_harmonise_interp(dar_tup):
    """ Perform linear interpolation on merged set of axis points for two or more xray DataArry objects.
    This function can be used to prepare (harmonise) multiple xray.DataArray objects for multiplication or addition
    on a common set of coordinate axis points by linearly interpolating all DataArray objects onto the same
    set of points, obtained by merging and sorting the points from all input DataArray objects.

    The DataArry objects provided.    The scipy linear grid interpolator is ued for this purpose. See:
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


