from mantid.simpleapi import logger
import numpy as np

__all__ = [
    "stitch_profiles",
]


def stitch_profiles(profiles, overlaps, target_profile_index=0):
    r"""
    Stitch together a sequence of intensity profiles with overlapping domains, returning a single encompassing profile.

    **drtsans objects used**:
    ~drtsans.dataobjects.IQmod
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/dataobjects.py>

    Parameters
    ----------
    profiles: list
        A list of  ~drtsans.dataobjects.IQmod objects, ordered with increasing Q-values
    overlaps: list
        A list of overlap regions in the shape (start_1, end_1, start_2, end_2, start_3, end_3,...).
    target_profile_index: int
        Index of the ``profiles`` list indicating the target profile, that is, the profile defining the final scaling.

    Returns
    -------
    ~drtsans.dataobjects.IQmod
    """
    # Guard clause to verify the profiles are ordered with increasing Q-values
    first_q_values = np.array([profile.mod_q[0] for profile in profiles])  # collect first Q-value for each profile
    if np.all(np.diff(first_q_values) > 0) is False:
        raise ValueError("The profiles are not ordered with increasing Q-values")

    # Guard clause to validate that the number of overlap boundaries is congruent with the number of intensity profiles
    if len(overlaps) != 2 * (len(profiles) - 1):
        raise ValueError("The number of overlaps is not appropriate to the number of intensity profiles")

    # Pair the overlaps into (start_q, end_q) pairs
    overlaps = [overlaps[i : i + 2] for i in range(0, len(overlaps), 2)]

    def scaling(target, to_target, starting_q, ending_q):
        r"""Utility function to find the scaling factor bringing the to_target profile to the target profile scaling"""
        # Find the data points of the "target" profile in the overlap region
        indexes_in_overlap = (target.mod_q > starting_q) & (target.mod_q < ending_q) & np.isfinite(target.intensity)
        q_values_in_overlap = target.mod_q[indexes_in_overlap]
        # Interpolate the "to_target" profile intensities at the previously found Q values
        good_values = np.isfinite(to_target.intensity)
        to_target_interpolated = np.interp(
            q_values_in_overlap,
            to_target.mod_q[good_values],
            to_target.intensity[good_values],
        )
        scale = sum(target_profile.intensity[indexes_in_overlap]) / sum(to_target_interpolated)
        if scale <= 0:
            raise ValueError(
                f"Scale number: {scale}. The scaling number for stitching cannot be negative. "
                + "Please check the stitching range or profile pattern"
            )
        else:
            return scale

    # We begin stitching to the target profile the neighboring profile with lower Q-values, then proceed until we
    # run out of profiles with lower Q-values than the target profile
    target_profile = profiles[target_profile_index]
    current_index = target_profile_index - 1
    while current_index >= 0:
        to_target_profile = profiles[current_index]
        start_q, end_q = overlaps[current_index]

        # Rescale the "to_target" profile to match the scaling of the target profile
        scale = scaling(target_profile, to_target_profile, start_q, end_q)
        to_target_profile = to_target_profile * scale
        print(f"Stitching profile {current_index} to profile {target_profile_index}. Scale factor is {scale:.3e}")

        # Discard extrema points
        to_target_profile = to_target_profile.extract(to_target_profile.mod_q < end_q)  # keep data with Q < end_q
        target_profile = target_profile.extract(target_profile.mod_q > start_q)  # keep data with Q > start_q

        # Stitch by concatenation followed by sorting, save the result into a new target profile
        target_profile = to_target_profile.concatenate(target_profile)  # just put one profile after the other
        target_profile = target_profile.sort()  # sort data points by increasing Q

        # Move to the next to-target profile
        current_index = current_index - 1

    # We continue stitching to the target profile the neighboring profile with higher Q-values, then proceed until we
    # run out of profiles with higher Q-values than the target profile
    current_index = target_profile_index + 1
    while current_index < len(profiles):
        to_target_profile = profiles[current_index]
        start_q, end_q = overlaps[current_index - 1]

        # Rescale the "to_target" profile to match the scaling of the target profile
        scale = scaling(target_profile, to_target_profile, start_q, end_q)
        to_target_profile = to_target_profile * scale
        print(f"Stitching profile {current_index} to profile {target_profile_index}. Scale factor is {scale:.3e}")

        # Discard extrema points
        to_target_profile = to_target_profile.extract(to_target_profile.mod_q > start_q)  # keep data with Q < end_q
        target_profile = target_profile.extract(target_profile.mod_q < end_q)  # keep data with Q > start_q

        # Stitch by concatenation followed by sorting, save the result into a new target profile
        target_profile = target_profile.concatenate(to_target_profile)  # just put one profile after the other
        target_profile = target_profile.sort()  # sort data points by increasing Q

        # Move to the next to-target profile
        current_index = current_index + 1

    return target_profile


def olt_q_boundary(reduction_config, iq1d_high_q, boundary, anisotropic=False):
    r"""Initialize the stitching boundaries when a list of boundaries has not been specified

    Parameters
    ----------
    reduction_config: dict
        The dictionary of reduction parameters
    iq1d_high_q: :py:obj:`~drtsans.dataobjects.IQmod`
         The unbinned higher-Q profile of the intensity profiles to stitch
    boundary: str
        Either 'min' or 'max'
    anisotropic: bool
        If True, returns the two boundaries for anisotropic/wedge binning, i.e. for "min" returns
        "wedge1overlapStitchQmin" and "wedge2overlapStitchQmin".
        If False, returns the boundary for isotropic binning, i.e. for "min" returns "overlapStitchQmin".

    Returns
    -------
    list
        A list of boundary values
    """
    if boundary not in ("min", "max"):
        raise ValueError('Only "min" or "max" are valid arguments')

    def get_isotropic_stitch_q_boundary():
        """Get stitch overlap boundary (min/max) for isotropic reduction

        Returns
        -------
        list
            List of one boundary value (min/max)
        """
        olt_q = reduction_config[f"overlapStitchQ{boundary}"]  # guaranteed `None` or `list`
        boundary_value = get_stitch_q_boundary(olt_q)
        return [boundary_value]

    def get_anisotropic_stitch_q_boundary():
        """Get stitch overlap boundary (min/max) for anisotropic reduction (wedges)

        Returns
        -------
        list
            List of one boundary value each (min/max) for the two wedges
        """
        boundary_values = []
        for wedge in ["wedge1", "wedge2"]:
            olt_q = reduction_config[f"{wedge}overlapStitchQ{boundary}"]  # guaranteed `None` or `list`
            boundary_value = get_stitch_q_boundary(olt_q)
            boundary_values.append(boundary_value)
        return boundary_values

    def get_stitch_q_boundary(overlap_q):
        """Get stitch overlap boundary (min/max) from configuration or boundary value of the higher-Q intensity profile

        Parameters
        ----------
        overlap_q: list or None
            The value of the configuration parameter for the overlap boundary

        Returns
        -------
        float
        """
        if overlap_q is None:
            logger.notice(f"Stitch Q{boundary} is None. Getting stitch Q{boundary} from {boundary}(I(Q))")
            extremum_function = getattr(iq1d_high_q.mod_q, boundary)  # either min() or max() method
            return extremum_function()
        else:
            # TODO: temporarily, until the midrange detector is added, the list only has one entry
            return overlap_q[0]

    if not anisotropic:  # scalar or annular
        return get_isotropic_stitch_q_boundary()
    else:  # wedges
        return get_anisotropic_stitch_q_boundary()
