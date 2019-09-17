import numpy as np
from drtsans.settings import namedtuplefy


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
    def find_n_truth(truths, repeats, reverse=False,
                     message='Could not find pixel index'):
        r"""
        Find first array index of consecutive `repeats` True values.

        Parameters
        ----------
        truths: list
            list of `True` and `False` items
        repeats: int
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
        truth_array = truths[::-1] if reverse else truths
        pat = [True]*repeats
        for i in range(len(truth_array) - repeats):
            if truth_array[i: i + repeats] == pat:
                return len(truths) - i - 1 if reverse else i
        else:
            raise IndexError(message)

    av = np.average(intensities)
    end_threshold = tube_threshold * av
    shadow_threshold = shadow_threshold * av

    # Find edges of the tube
    illuminated = [bool(i > end_threshold) for i in intensities]
    bottom_pixel = find_n_truth(illuminated, tube_edge_min_width,
                                message='Could not find bottom tube edge')
    top_pixel = find_n_truth(illuminated, tube_edge_min_width,
                             message='Could not find top tube edge',
                             reverse=True)

    # Find the shadow region
    active = intensities[bottom_pixel: top_pixel - tube_edge_min_width]
    shadowed = [bool(i < shadow_threshold) for i in active]
    bottom_shadow_pixel = bottom_pixel +\
        find_n_truth(shadowed, shadow_edge_min_width,
                     message='Could not find bottom shadow edge')

    active = intensities[bottom_shadow_pixel + shadow_edge_min_width:
                         top_pixel - tube_edge_min_width]
    illuminated = [bool(i > shadow_threshold) for i in active]
    above_shadow_pixel = bottom_shadow_pixel + shadow_edge_min_width +\
        find_n_truth(illuminated, 1,
                     message='could not find above shadow region')

    # Check for a faulty tube
    active_tube_length = top_pixel - bottom_pixel + 1
    shadow_length = above_shadow_pixel - bottom_shadow_pixel
    if active_tube_length < min_illuminated_length + shadow_length:
        raise RuntimeError('Faulty tube found')

    return dict(bottom_pixel=bottom_pixel, top_pixel=top_pixel,
                bottom_shadow_pixel=bottom_shadow_pixel,
                above_shadow_pixel=above_shadow_pixel)
