# local imports
from drtsans.dataobjects import IQmod, I1DAnnular

# third party imports
from mantid.simpleapi import logger
import numpy as np

# standard imports


__all__ = [
    "stitch_profiles",
]


def stitch_profiles(profiles, overlaps, target_profile_index=0):
    r"""
    Stitch together a sequence of intensity profiles with overlapping domains, returning a single encompassing profile.

    **drtsans objects used**:
    ~drtsans.dataobjects.IQmod
    <https://github.com/neutrons/drtsans/blob/next/src/drtsans/dataobjects.py>

    Parameters
    ----------
    profiles: list
        A list of  ~drtsans.dataobjects.IQmod objects, ordered with increasing Q-values
    overlaps: list of lists or list
        The overlap regions either as: [(start_1, end_1), (start_2, end_2), (start_3, end_3), ...] or (for backwards
        compatibility): [start_1, end_1, start_2, end_2, start_3, end_3, ...]
    target_profile_index: int
        Index of the ``profiles`` list indicating the target profile, that is, the profile defining the final scaling.

    Returns
    -------
    ~drtsans.dataobjects.IQmod

    Raises
    -------
    ValueError
        If either the arguments are incorrect ((i) profiles not in order or increasing Q or (ii) the number of overlaps
        not congruent with the number of profiles or (iii) overlaps is not a list of numbers or list of
        lists/tuples) or a stitching scaling factor <= 0 is calculated.
    """
    # Guard clause to verify the profiles are ordered with increasing Q-values
    first_q_values = np.array([profile.mod_q[0] for profile in profiles])  # collect first Q-value for each profile
    if np.all(np.diff(first_q_values) > 0) is False:
        raise ValueError("The profiles are not ordered with increasing Q-values")

    # Check overlaps parameter and if needed convert from list of numbers to list of pairs of numbers
    if not overlaps:
        raise ValueError("Must provide list of overlaps")
    if all(isinstance(x, (int, float)) for x in overlaps) and len(overlaps) == 2 * (len(profiles) - 1):
        # Convert to new format by pairing the overlaps into (start_q, end_q) pairs
        overlaps = [overlaps[i : i + 2] for i in range(0, len(overlaps), 2)]
    if not all(isinstance(x, (list, tuple)) and len(x) == 2 for x in overlaps):
        raise ValueError("Parameter overlaps must be either list of numbers or list of lists of numbers")

    # Guard clause to validate that the number of overlap boundaries is congruent with the number of intensity profiles
    if len(overlaps) != len(profiles) - 1:
        raise ValueError("The number of overlaps is not appropriate to the number of intensity profiles")

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


def get_stitch_boundaries(reduction_config, iq1d_unbinned, anisotropic=False):
    r"""Initialize the stitching boundaries when a list of boundaries has not been specified

    Parameters
    ----------
    reduction_config: dict
        The dictionary of reduction parameters
    iq1d_unbinned: list of :py:obj:`~drtsans.dataobjects.IQmod`
         The unbinned IQmod profiles, one for each detector.
    anisotropic: bool
        If True, two binnings exist and the applicable boundaries are "wedge1overlapStitchMin",
        "wedge1overlapStitchMax", "wedge2overlapStitchMin", and "wedge2overlapStitchMax".
        If False, one binning exists and the applicable boundaries are "overlapStitchMin" and "overlapStitchMax".

    Returns
    -------
    list
        A list of lists of boundary values, for example for anisotropic=True:
        [
            [(start_overlap1, end_overlap1), (start_overlap2, end_overlap2)],  # wedge 1
            [(start_overlap1, end_overlap1), (start_overlap2, end_overlap2)]   # wedge 2
        ]
        For anisotropic=False:
        [
            [(start_overlap1, end_overlap1), (start_overlap2, end_overlap2)]
        ]
    """
    bounds = []

    def profile_boundaries(qmin_values, qmax_values):
        r"""

        Parameters
        ----------
        qmin_values: None or list
            A list of Qmin values to use, one per stitching overlap region. For py:obj:`None`, uses the min value of
            the higher-Q profile of each pair of profiles to stitch.
        qmax_values: None or list
            A list of Qmax values to use, one per stitching overlap region. For py:obj:`None`, uses the max value of
            the higher-Q profile of each pair of profiles to stitch.

        Returns
        -------
        list
            A list of boundary values [(start, end), (start,end), ...]
        """
        if not qmin_values or not qmax_values:  # catches None or empty list
            # use min/max from the Q range for the higher Q profile, e.g. for [p1, p2, p3]
            # stitching boundaries for [p1, p2] will come from min(p2), max(p2)
            # stitching boundaries for [p2, p3] will come from min(p3), max(p3)
            qmin_values = [x.mod_q.min() for x in iq1d_unbinned[1:]]
            qmax_values = [x.mod_q.max() for x in iq1d_unbinned[1:]]
            logger.notice("Stitch Qmin/max is None. Getting stitch boundaries from min(I(Q)) and max(I(Q))")
        profile_bounds = []
        for qmin, qmax in zip(qmin_values, qmax_values):
            profile_bounds.append((qmin, qmax))
        return profile_bounds

    if not anisotropic:
        qmin_values = reduction_config["overlapStitchQmin"]
        qmax_values = reduction_config["overlapStitchQmax"]
        bounds.append(profile_boundaries(qmin_values, qmax_values))
    else:
        for wedge in ["wedge1", "wedge2"]:
            qmin_values = reduction_config[f"{wedge}overlapStitchQmin"]
            qmax_values = reduction_config[f"{wedge}overlapStitchQmax"]
            bounds.append(profile_boundaries(qmin_values, qmax_values))

    return bounds


def get_target_profile_index(iq1d_binned, reduction_config):
    """Get the target profile index using the binned IQ1D data and reduction configuration.

    The reduction configuration field ``overlapStitchReferenceDetector`` defines the name of the detector to use.
    The size of ``iq1d_binned`` determines the ordering of the data, and hence can be used to return the index.
    If ``len(iq1d_binned)`` == 2, its order is [main, wing]
    If ``len(iq1d_binned)`` == 3, its order is [main, midrange, wing]

    Parameters
    ----------
    iq1d_binned: list of lists of ~drtsans.dataobjects.IQmod
        A list of lists of ~drtsans.dataobjects.IQmod objects. The outer list is the list of detectors ordered by
        increasing Q-values. The inner lists are for different binnings of the original data from the detector.
        Example: [[IQ_main_wedge1, IQ_main_wedge2], [IQ_wing_wedge1, IQ_wing_wedge2]]
    reduction_config: dict
        The dictionary of reduction parameters.

    Returns
    -------
    int
        The index of the target profile
    """
    detector_index_map = {
        2: {
            "main": 0,
            "wing": 1,
        },
        3: {
            "main": 0,
            "midrange": 1,
            "wing": 2,
        },
    }

    detector_key = reduction_config["overlapStitchReferenceDetector"]
    return detector_index_map[len(iq1d_binned)][detector_key]


def stitch_binned_profiles(iq1d_unbinned, iq1d_binned, reduction_config):
    """Stitch together sequences of intensity profiles from different detector panels (main, wing and midrange),
    that cover a different range of Q-values, returning a single encompassing profile per sequence.

    When "1DQbinType" is "scalar" or "annular", there is one sequence of intensity profiles, i.e. one intensity profile
    binning per detector panel, and one combined intensity profile is returned.
    When "1DQbinType" is "wedge", there are two sequences of intensity profiles, i.e. two intensity profile binnings
    per detector panel, and two combined intensity profiles are returned.

    Parameters
    ----------
    iq1d_unbinned: list
        A list of ~drtsans.dataobjects.IQmod objects for different detectors. The intensity profiles to stitch together
        ordered by increasing Q-values.
    iq1d_binned: list of lists of ~drtsans.dataobjects.IQmod
        A list of lists of ~drtsans.dataobjects.IQmod objects. The outer list is the list of detectors ordered by
        increasing Q-values. The inner lists are for different binnings of the original data from the detector.
        Example: [[IQ_main_wedge1, IQ_main_wedge2], [IQ_wing_wedge1, IQ_wing_wedge2]]
    reduction_config: dict
        The dictionary of reduction parameters.

    Returns
    -------
    list of ~drtsans.dataobjects.IQmod objects
        The list of combined (stitched) profiles. When "1DQbinType" is "scalar" or "annular", the list contains one
        combined intensity profile. When "1DQbinType" is "wedge", the list contains two combined intensity profiles,
        one per wedge.
    """
    iq1d_combined_out = []

    if isinstance(iq1d_binned[0][0], I1DAnnular):
        # stitching of annular profiles is not supported, return empty combined profile(s)
        iq1d_binned_main = iq1d_binned[0]
        for ibinning in range(len(iq1d_binned_main)):
            iq1d_combined = I1DAnnular(intensity=[], error=[], phi=[])
            iq1d_combined_out.append(iq1d_combined)
        logger.notice("Skipping stitching. Stitching of annular profiles is not supported.")
        return iq1d_combined_out

    boundaries = get_stitch_boundaries(reduction_config, iq1d_unbinned, anisotropic=len(iq1d_binned[0]) > 1)

    iq1d_binned_main = iq1d_binned[0]
    for ibinning in range(len(iq1d_binned_main)):
        profiles = [detector_profiles[ibinning] for detector_profiles in iq1d_binned]
        overlaps = boundaries[ibinning]
        # do stitching
        try:
            iq1d_combined = stitch_profiles(
                profiles=profiles,
                overlaps=overlaps,
                target_profile_index=get_target_profile_index(iq1d_binned, reduction_config),
            )
        except ValueError:
            iq1d_combined = IQmod(intensity=[], error=[], mod_q=[], delta_mod_q=[])
        iq1d_combined_out.append(iq1d_combined)
    return iq1d_combined_out
