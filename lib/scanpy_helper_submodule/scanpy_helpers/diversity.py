import numpy as np
from numba import njit


@njit
def gini_index(array: np.ndarray) -> float:
    """
    Calculate the Gini coefficient of a numpy array.

    TODO fix docstring and doctests

    Based on: https://github.com/oliviaguest/gini
    Args:
        array (array-like): input array
    Returns:
        float: gini-index of ``array``

    >>> a = np.zeros((10000))
    >>> a[0] = 1.0
    >>> '%.3f' % gini(a)
    '1.000'
    >>> a = np.ones(100)
    >>> '%.3f' % gini(a)
    '0.000'
    >>> a = np.random.uniform(-1,0,1000000)
    >>> '%.2f' % gini(a)
    '0.33'
    """
    if not np.all(array >= 0):
        raise ValueError("Gini index is only defined for non-negative values. ")
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array += 1e-12  # values cannot be 0
    array = np.sort(array)  # values must be sorted
    index = np.arange(1, array.shape[0] + 1)  # index per array element
    n = array.shape[0]  # number of array elements
    return (np.sum((2 * index - n - 1) * array)) / (
        n * np.sum(array)
    )  # Gini coefficient
