import numpy as np

r""" Hyperlinks to mantid algorithms
CloneWorkspace <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
DeleteWorkspaces <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspaces-v1.html>
Integration <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
MaskDetectors <https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html>
MaskDetectorsIf <https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html>
ReplaceSpecialValues <https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html>
"""
from mantid.simpleapi import (CloneWorkspace, DeleteWorkspaces, Integration, MaskDetectors, MaskDetectorsIf,
                              ReplaceSpecialValues)
from mantid.api import mtd

r"""
Hyperlinks to drtsans functions
namedtuplefy, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.tubecollection import TubeCollection


__all__ = ['apparent_tube_width', 'find_edges', 'fit_positions']


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


def apparent_tube_width(input_workspace, output_workspace=None):
    r"""
    Determine the tube widhth most efficient for detecting neutrons. An effective tube (or pixel) diameter is
    determined for tubes in the front panel, and likewise for the tubes in the back panel.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Input workspace, usually a flood run.
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of ``input_workspace`` is used, thus
        calibrating the pixel widths of the input workspace.

    **Mantid algorithms used:**
        :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,
        :ref:`DeleteWorkspaces <algm-DeleteWorkspaces-v1>`,
        :ref:`Integration <algm-Integration-v1>`,
        :ref:`MaskDetectors <algm-MaskDetectors-v1>`,
        :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
        :ref:`ReplaceSpecialValues <algm-ReplaceSpecialValues-v1>`,

    Returns
    -------
    ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    integrated_intensities = unique_workspace_dundername()
    Integration(InputWorkspace=input_workspace, OutputWorkspace=integrated_intensities)

    # Mask non-finite intensities (nan, inf). They can't be used in the calculation.
    #
    # Replace non-finite intensities with a value of -1
    ReplaceSpecialValues(InputWorkspace=integrated_intensities, OutputWorkspace=integrated_intensities,
                         NanValue=-1, NanError=-1, InfinityValue=-1, InfinityError=-1)
    # Mask detectors with negative intensities
    mask_workspace = unique_workspace_dundername()
    MaskDetectorsIf(InputWorkspace=integrated_intensities, Operator='Less', Value=0., OutputWorkspace=mask_workspace)
    MaskDetectors(Workspace=integrated_intensities, MaskedWorkspace=mask_workspace)

    # Calculate the count density for each tube. Notice that if the whole tube is masked, then the associated
    # intensity is stored as nan.
    #
    # Sort the tubes according to the X-coordinate in decreasing value. This is the order when sitting on the
    # sample and iterating over the tubes "from left to right"
    collection = TubeCollection(integrated_intensities, 'detector1').sorted(view='decreasing X')
    intensities = mtd[integrated_intensities].extractY().flatten()
    count_densities = list()
    for tube in collection:
        weighted_intensities = tube.readY.ravel() / tube.pixel_heights
        d = np.mean(weighted_intensities[~tube.isMasked])
        count_densities.append(d)
    count_densities = np.array(count_densities)  # is convenient to cast densities into a numpy array data structure.

    # Determine the count densities per panel and for the whole detector array.
    # We must be careful to pick only tubes with finite densities (avoid 'nan')
    average_count_density = np.mean(count_densities[np.isfinite(count_densities)])
    front_count_density = np.mean(count_densities[::2][np.isfinite(count_densities[::2])])  # front tubes, even indexes
    back_count_density = np.mean(count_densities[1::2][np.isfinite(count_densities[1::2])])  # back tubes, odd indexes

    # Determine the front and back pixel widths, then compare to test data
    nominal_width = collection[0][0].width  # width of the first pixel in the first tube
    front_width = (front_count_density / average_count_density) * nominal_width
    back_width = (back_count_density / average_count_density) * nominal_width

    # Insert the updated pixel widths in the output workspace.
    if output_workspace != str(input_workspace):  # are we overwriting the pixel widths of the input workspace?
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)
    collection = TubeCollection(output_workspace, 'detector1').sorted(view='decreasing X')
    for tube in collection[::2]:  # front tubes
        for pixel in tube:
            pixel.width = front_width
    for tube in collection[1::2]:  # back tubes
        for pixel in tube:
            pixel.width = back_width

    DeleteWorkspaces(integrated_intensities, mask_workspace)
    return mtd[output_workspace]
