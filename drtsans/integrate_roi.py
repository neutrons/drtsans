from mantid.api import AnalysisDataService as mtd


def integrate_detector_roi(workspace, roi_det_list):
    """ Integrate neutron counts in ROI of detectors (pixels)
    Parameters
    ----------
    workspace: String or ~mantid.api.MatrixWorkspace
        Name of workspace or Workspace instance
    roi_det_list: List of integers
        Workspace indexes for the detectors in ROI

    Returns
    -------
    Float
        Integrated intensity
    """
    # Get workspace instance if input is a string
    workspace = mtd[str(workspace)]

    # Get Y array
    counts_array = workspace.extractY()

    # Get the sum
    roi_intensity = counts_array[roi_det_list]
    print('DDD:', roi_intensity)
    roi_intensity = roi_intensity.sum()

    return roi_intensity
