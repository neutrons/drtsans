import numpy as np

# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py
from drtsans.settings import namedtuplefy


def _consecutive_true_values(values, how_many, reverse=False,
                             message='Could not find pixel index'):
    r"""
    Find first array index of consecutive `how_many` True values.

    Parameters
    ----------
    values: list
        list of `True` and `False` items
    how_many: int
        Number of desired consecutive `True` values
    message: str
        Exception message

    Returns
    -------
    int

    Raises
    ------
    IndexError
        If no index is found for any of the edges
    RuntimeError
        If a faulty tube is found
    """
    # use the array or the reverse one
    truth_array = values[::-1] if reverse else values
    # create a sub-array of length how_many of True values that we want to find
    pattern = [True]*how_many
    # loop over the input data and return the first index where the next
    # how_many elements match the pattern
    for i in range(len(truth_array) - how_many):
        if truth_array[i: i + how_many] == pattern:
            return len(values) - i - 1 if reverse else i
    # raise an error if the pattern is not found
    else:
        raise IndexError(message)


@namedtuplefy
def find_edges(intensities, tube_threshold=0.2, shadow_threshold=0.3,
               tube_edge_min_width=3, shadow_edge_min_width=4,
               min_illuminated_length=7):
    r"""
    Find the active length of the tube and the shadow region

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    Parameters
    ----------
    intensities: list
        pixel intensities along the tube.
    tube_threshold: float
        fraction of the average intensity to determine the tube edges.
    shadow_threshold: float
        fraction of the average intensity to determine the shadow region.
    tube_edge_min_width: int
        required minimum number of consecutive good pixels above
        the tube threshold
    shadow_edge_min_width: int
        required minimum number of consecutive shadowed pixels
    min_illuminated_length: int
        minimum number of illuminated pixels on the active length

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - bottom_pixel: first illuminated pixel
        - top_pixel: last illuminated pixel
        - bottom_shadow_pixel: first shadowed pixel
        - above_shadow_pixel= first illuminated pixel above the shadow region
    """
    # calculate minimum intensity thresholds for tube ends and shaddows
    av = np.average(intensities)
    end_threshold = tube_threshold * av
    shadow_threshold = shadow_threshold * av

    # Find edges of the tube: want at least tube_edge_min_width pixels
    # (starting from the top or bottom of a tube) that have intensities greater
    # than the threshold
    illuminated = [bool(i > end_threshold) for i in intensities]
    bottom_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                            message='Could not find bottom tube edge')
    top_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                         message='Could not find top tube edge',
                                         reverse=True)

    # Find the shadow region: similar to tube edges, but in this case
    # we want shadow_edge_min_width intensities less than the shadow threshold,
    # followed by at least one intensity greater than the threshold
    active = intensities[bottom_pixel: top_pixel - tube_edge_min_width]
    shadowed = [bool(i < shadow_threshold) for i in active]
    bottom_shadow_pixel = bottom_pixel +\
        _consecutive_true_values(shadowed, shadow_edge_min_width,
                                 message='Could not find bottom shadow edge')

    active = intensities[bottom_shadow_pixel + shadow_edge_min_width:
                         top_pixel - tube_edge_min_width]
    illuminated = [bool(i > shadow_threshold) for i in active]
    above_shadow_pixel = bottom_shadow_pixel + shadow_edge_min_width +\
        _consecutive_true_values(illuminated, 1,
                                 message='could not find above shadow region')

    # Check for a faulty tube: we want a certain number of pixels not in the bar shaddow
    active_tube_length = top_pixel - bottom_pixel + 1
    shadow_length = above_shadow_pixel - bottom_shadow_pixel
    if active_tube_length < min_illuminated_length + shadow_length:
        raise RuntimeError('Faulty tube found')

    return dict(bottom_pixel=bottom_pixel, top_pixel=top_pixel,
                bottom_shadow_pixel=bottom_shadow_pixel,
                above_shadow_pixel=above_shadow_pixel)


@namedtuplefy
def fit_positions(edge_pixels, bar_positions, tube_pixels=256):
    r"""
    Fit the position and heights of the pixels in a tube. The bar_positions as a function of
    edge pixels are fitted to a 5th degree polynomial. The positions of the pixels along the
    tube are the values of the polynomial at integer points, while the heights are the derivatives.

    Description from the master requirements document, section A2.1

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    Uses :ref:`~numpy.polynomial.polynomial.polyfit`.

    Parameters
    ----------
    edge_pixels: list (or numpy array)
        the bottom pixel for each bar position, as found in `find_edges` function
    bar_positions: list (or numpy array)
        the bar position from the logs for each file in the bar scan
    tube_pixels: integer
        number of pixels for which to calculate positions and heights

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - calculated_positions: calculated positions of the pixels
        - calculated_heights: calculated pixel heights
    """
    message_len = 'The positions of the bar and edge pixels have to be the same length'
    assert len(edge_pixels == bar_positions), message_len

    # fit the bar positions to a 5th degree polynomial in edge_pixels
    coefficients = np.polynomial.polynomial.polyfit(edge_pixels, bar_positions, 5)
    # calculate the coefficients of the derivative
    deriv_coefficients = np.polynomial.polynomial.polyder(coefficients)
    # evalutae the positions
    calculated_positions = np.polynomial.polynomial.polyval(np.arange(tube_pixels), coefficients)
    # evaluate the heights
    calculated_heights = np.polynomial.polynomial.polyval(np.arange(tube_pixels), deriv_coefficients)

    return dict(calculated_positions=calculated_positions, calculated_heights=calculated_heights)
