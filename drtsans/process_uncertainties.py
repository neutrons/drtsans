from drtsans.samplelogs import SampleLogs
from mantid.simpleapi import mtd, SetUncertainties, MaskBins
import numpy


def set_init_uncertainties(input_workspace, output_workspace=None, mask_band_gap=True):
    """
    Set the initial uncertainty of a :py:obj:`~mantid.api.MatrixWorkspace`

    Mantid algorithm :ref:`SetUncertainties <algm-SetUncertainties-v1>` will be called to make sure
    1: set the uncertainty to square root of intensity
    2: make sure all zero uncertainties will be set to 1

    In case of output workspace is py:obj:`None`, the input workspace will be
    replaced by output workspace.

    If the workspace is in unit of Wavelength, all intensities and uncertainties within band gap shall be masked

    :exception RuntimeError: output workspace (string) is empty

    **Mantid algorithms used:**
    :ref:`SetUncertainties <algm-SetUncertainties-v1>`

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace
        Input workspace
    output_workspace: str
        Output workspace (workspace name or instance) or py:obj:`None` for in-place operation
    mask_band_gap: bool
        True if workspace is in Wavelength unit

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    input_workspace = str(input_workspace)
    if output_workspace is None:
        output_workspace = input_workspace  # in-place by default
    else:
        output_workspace = str(output_workspace)

    # Calculate uncertainties as square root and set 1 for 0 count
    # But SetUncertainties does not treat nan as SANS team desires
    SetUncertainties(InputWorkspace=input_workspace,
                     OutputWorkspace=output_workspace,
                     SetError='sqrtOrOne')

    # get a handle to the workspace
    output_ws = mtd[output_workspace]

    # Set nan as the uncertainty for all nan-intensity
    for ws_index in range(output_ws.getNumberHistograms()):
        vec_y = output_ws.readY(ws_index)
        nan_indexes = numpy.argwhere(numpy.isnan(vec_y))

        # There existing nan
        if len(nan_indexes) > 0:
            vec_e = output_ws.dataE(ws_index)
            vec_e[nan_indexes] = numpy.nan
        # END-IF
    # END-FOR (spectra)

    # Set band gap
    if mask_band_gap and output_ws.getAxis(0).getUnit().unitID() == 'Wavelength':
        output_ws = _mask_bins_in_band_gap(output_ws)

    return output_ws


def _mask_bins_in_band_gap(workspace):
    """
    If data is measured at 30Hz and logs lead_max and skip_min
    :exception  RuntimeError: if lead_max or skip_min cannot be found

    **Mantid algorithms used:**
    :ref:`MaskBins <algm-MaskBins-v1>`

    Parameters
    ----------
    workspace: ~mantid.api.MatrixWorkspace
        Workspace to correct. Always done inplace

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    # use sample log object to simplify this
    samplelogs = SampleLogs(workspace)

    # Use BL6:Chop:Skf1:SpeedUserReq to get chopper frequency to check only applied to chopper @ 30 Hz
    chopper_frequency = samplelogs['BL6:Chop:Skf1:SpeedUserReq'].value.mean()
    if abs(chopper_frequency - 30) > 5.:
        # if chopper is not 30 Hz
        return workspace

    # Check sample logs
    lead_max_wl = samplelogs['lead_max'].value
    skip_min_wl = samplelogs['skip_min'].value

    # Mask a range of X-values
    workspace = MaskBins(InputWorkspace=workspace, OutputWorkspace=workspace, XMin=lead_max_wl, XMax=skip_min_wl)

    return workspace
