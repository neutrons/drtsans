import numpy as np


def stitch_profiles(profiles, overlaps, target_profile=0):
    r"""
    Stitch together a sequence of intensity profiles with overlapping domains, returning a single encompassing profile.

    Parameters
    ----------
    profiles: list
        A list of  ~drtsans.dataobjects.IQmod objects.
    overlaps: list
        A list of overlap regions in the shape (start_1, end_1, start_2, end_2, start_3, end_3,...)
    target_profile: int
        Index of the ``profiles`` list indicating which profile sets the final scaling.

    Returns
    -------
    ~drtsans.dataobjects.IQmod
    """
    # Guard clause to validate that the number of overlap boundaries is congruent with the number of intensity profiles
    if len(overlaps) != 2 * (len(profiles) - 1):
        raise ValueError('The number of overlaps is not appropriate to the number of intensity profiles')

    final_scaling = 1  # scaling that must be used to scale the stitched profile to the scale of the desired profile
    low_q_profile = profiles[0]
    # low_q_profile is an IQmod object with the following attributes:
    #   low_q_profile.intensity: array of intensities
    #   low_q_profile.error: array of uncertainties for the intensities
    #   low_q_profile.mod_q: array of Q values
    #   low_q_profile.delta_mod_q: array of uncertainties for the Q values
    for i in range(1, len(profiles)):  # iterate, stitching two profiles at each iteration
        start_q, end_q = overlaps.pop(0), overlaps.pop(0)  # pick the first overal region from the list of overlaps
        high_q_profile = profiles[i]  # another IQmod object

        # Find the data points of the low-Q profile in the overlap region
        indexes_in_overlap = (low_q_profile.mod_q > start_q) & (low_q_profile.mod_q < end_q)
        q_values_in_overlap = low_q_profile.mod_q[indexes_in_overlap]

        # Interpolate the high-Q profile intensities at the previously found Q values
        high_q_interpolated = np.interp(q_values_in_overlap, high_q_profile.mod_q, high_q_profile.intensity)

        # Rescale the high-Q profile to match the scaling of the low-Q profile
        scaling = sum(low_q_profile.intensity[indexes_in_overlap]) / sum(high_q_interpolated)
        high_q_profile = scaling * high_q_profile

        # Discard extrema points
        low_q_profile = low_q_profile.select(low_q_profile.mod_q < end_q)  # keep data with Q < end_q
        high_q_profile = high_q_profile.select(high_q_profile.mod_q > start_q)  # keep data with Q > start_q

        # Stitch by concatenation followed by sorting, save the result into a new low_q_profile
        low_q_profile = low_q_profile.concatenate(high_q_profile)  # just put one profile after the other
        low_q_profile = low_q_profile.sort()  # sort data points by increasing Q

        if i == target_profile:  # we found the profile to which we want to scale all others
            final_scaling = 1. / scaling

    return final_scaling * low_q_profile
