"""Utility functions.

By: Genomic Data Modeling Lab
"""

import numbers
from importlib import metadata
import datetime
import os
import re
import numpy as np

NUMPY_NUMERIC_DTYPE_KINDS = ("f", "u", "i")


def is_numeric_nparray(x):
    """Test whether numpy array is numeric.

    Note: that a mixture of booleans and numeric may be

    """
    return x.dtype.kind in NUMPY_NUMERIC_DTYPE_KINDS


# TODO what to do about np.nan
# TODO is it OK that object arrays raise TypeError as opposed to
# returning false, see unit test for example

def is_biallelic(x):
    """Test whether data set is biallelic (0,1, np.nan) numeric."""
    data_set = np.unique(x, equal_nan=True)

    if data_set.size > 3 or not is_numeric_nparray(data_set):
        return False

    delta = np.setdiff1d(data_set, np.array([0,1]))

    return delta.size == 0 or np.isnan(delta).all()


def get_version(file=None):
    return metadata.version("afcn")


def is_int(string_value):
    """Check whether input string is an integer."""
    return re.match("^[+-]?\d+$", string_value) is not None


def is_float(string_value):
    """Check whether input string is a floating point number."""
    return re.match("^[+-]?\d+\.\d*$", string_value) is not None
