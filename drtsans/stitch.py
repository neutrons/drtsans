import collections
import numpy as np

r"""
Hyperlinks to Mantid algorithms
CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
DeleteWorkspace <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v3.html>
Stitch1D <https://docs.mantidproject.org/nightly/algorithms/Stitch1D-v3.html>
"""
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, Stitch1D

r"""
Hyperlinks to drtsans functions
unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
IQmod <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/dataobjects.py>
"""  # noqa: E501
from drtsans.settings import unique_workspace_dundername
from drtsans.dataobjects import IQmod


def stitch_intensities(intensities=None, overlaps=None):
    r"""
    Stitch together a sequence of intensity profiles with overlapping domains to return a single encompassing profile.

    Parameters
    ----------
    intensities: list
        A list of I(Q) intensities. You can pass objets of type ~drtsans.dataobjects.IQmod and/or all the four
        components separately and in this order: intensity, error, mod_q, delta_mod_q. Each of these components must be
        a one-dimensional array.
        Allowed syntax for this list:
        - (intensity_1, error_1 mod_q_1 delta_mod_q_1, intensity_2 error_2 mod_q_2 delta_mod_q_2,...)
        - (IQmod_1, IQmod_2, IQmod_3,...)
        - (IQmod_1, intensity_2 error_2 mod_q_2 delta_mod_q_2, IQmod_3,...)
    overlaps: list
        A list of overlap regions. Allowed syntax for this list:
        - ((start_1, end_1), (start_2, end_2), (start_3, end_3),...)
        - (start_1, end_1, start_2, end_2, start_3, end_3,...)

    **Mantid algorithms used:**
    :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`,
    :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`,
    :ref:`Stitch1D <algm-Stitch1D-v3>`,

    Returns
    -------
    ~drtsans.dataobjects.IQmod
    """
    # We convert the list of intensities to a list of IQmod objects, i.e [IQmod_1, IQmod_2,...]
    # We need to iterate the list of intensities and whenever we found an item that it's not an IQmod object, we
    # create one by replacing the four components of the intensity profile.
    intensities = list(intensities)  # convert to list, necessary if passing a tuple
    index = 0
    while len(intensities) > index:
        if isinstance(intensities[index], IQmod) is False:
            # we found the four components: intensity, error, mod_q, and delta_mod_q
            new_intensity_profile = IQmod(*intensities[index: index+4])
            intensities = intensities[0:index] + [new_intensity_profile] + intensities[index+4:]
        index += 1

    # We convert the list of overlap regions to a single list of alternating boundaries,
    # i.e. [start_1, end_1, start_2, end_2,...]
    overlaps = list(overlaps)  # convert to list, necessary if passing a tuple
    if isinstance(overlaps[0], collections.Sequence):  # we're passing ((start_1, end_1), (start_2, end_2),...)
        overlaps = [boundary for overlap_region in overlaps for boundary in overlap_region]

    # Guard clause to validate that the number of overlap boundaries is congruent with the number of intensity profiles
    if len(overlaps) != 2 * (len(intensities) - 1):
        raise ValueError('The number of overlaps is not appropriate to the number of intensity profiles')

    # Stitch together the first two intensity profiles
    i_of_q_1, i_of_q_2 = intensities.pop(0), intensities.pop(0)
    start_overlap, end_overlap = overlaps.pop(0), overlaps.pop(0)

    # first interpolate the data in the overlap region
    interpolation_domain = np.linspace(start_overlap, end_overlap, 100)  # will 100 be always enough?
    low_q_workspace = CreateWorkspace(DataX=interpolation_domain, NSpec=1, EnableLogging=False,
                                      DataY=np.interp(interpolation_domain, i_of_q_1.mod_q, i_of_q_1.intensity),
                                      OutputWorkspace=unique_workspace_dundername())
    high_q_workspace = CreateWorkspace(DataX=interpolation_domain, NSpec=1, EnableLogging=False,
                                       DataY=np.interp(interpolation_domain, i_of_q_2.mod_q, i_of_q_2.intensity),
                                       OutputWorkspace=unique_workspace_dundername())

    # Stitch with Mantid's Stich1D to find out by how much we have to scale the high-Q profile
    _, scale = Stitch1D(LHSWorkspace=low_q_workspace, RHSWorkspace=high_q_workspace,
                        StartOverlap=start_overlap, EndOverlap=end_overlap)

    # Fuse the two intensity profiles as the concatenation of low_q_spectrum[q<start_overlap]
    # and high_q_spectrum[q>start_overlap]
    low_q_range = i_of_q_1.mod_q < start_overlap
    high_q_range = i_of_q_2.mod_q > start_overlap
    merged_q_range = np.concatenate((i_of_q_1.mod_q[low_q_range], i_of_q_2.mod_q[high_q_range]))
    stitched_intensity = np.concatenate((i_of_q_1.intensity[low_q_range], scale*i_of_q_2.intensity[high_q_range]))
    error_in_intensity = np.concatenate((i_of_q_1.error[low_q_range], scale*i_of_q_2.error[high_q_range]))
    error_in_q = np.concatenate((i_of_q_1.delta_mod_q[low_q_range], i_of_q_2.delta_mod_q[high_q_range]))
    new_intensity_profile = IQmod(stitched_intensity, error_in_intensity, merged_q_range, error_in_q)

    # clean up temporary workspaces
    [workspace.delete() for workspace in (low_q_workspace, high_q_workspace)]

    if len(intensities) > 0:  # there are intensity profiles left that need stitching
        # insert the stitched profile in the list of intensities at the top of the list
        intensities.insert(0, IQmod(stitched_intensity, error_in_intensity, merged_q_range, error_in_q))
        # A recursive call will stitch the newly create stitched profile to the next intensity profile down the]
        # list of intensities. Note that at this moment we have reduced the number of intensity profiles by one.
        stitch_intensities(intensities, overlaps)
    else:
        return new_intensity_profile
