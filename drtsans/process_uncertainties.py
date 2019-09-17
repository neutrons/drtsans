from mantid.simpleapi import SetUncertainties, MaskBins
import numpy


def set_init_uncertainties(input_ws, output_ws=None, mask_band_gap=True):
    """
    Set the initial uncertainty of a MatrixWorkspace
    Mantid algorithm SetUncertainties will be called to make sure
    1: set the uncertainty to square root of intensity
    2: make sure all zero uncertainties will be set to 1

    In case of output workspace is None, the input workspace will be
    replaced by output workspace.

    If the workspace is in unit of Wavelength, all intensities and uncertainties within band gap shall be masked

    :exception RuntimeError: output workspace (string) is empty

    :param input_ws: Input workspace
    :param output_ws: Output workspace (workspace name or instance) or None
    :param mask_band_gap: True if workspace is in Wavelength unit
    :return: reference to output workspace
    """
    # Set output workspace
    if output_ws is None:
        output_ws_name = str(input_ws)
        in_out_same = True
    else:
        output_ws_name = str(output_ws).strip()
        in_out_same = False
        if output_ws_name == '':
            raise RuntimeError('Output workspace name cannot be an empty '
                               'string')

    # Calculate uncertainties as square root and set 1 for 0 count
    # But SetUncertainties does not treat nan as SANS team desires
    output_ws = SetUncertainties(InputWorkspace=input_ws,
                                 OutputWorkspace=output_ws_name,
                                 SetError='sqrtOrOne')

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
    if in_out_same:
        input_ws = output_ws
    if mask_band_gap and input_ws.getAxis(0).getUnit().unitID() == 'Wavelength':
        output_ws = _mask_bins_in_band_gap(input_ws, output_ws)

    return output_ws


def _mask_bins_in_band_gap(original_ws, output_ws):
    """
    If data is measured at 30Hz and logs lead_max and skip_min
    :exception  RuntimeError: if lead_max or skip_min cannot be found
    :param original_ws:
    :param output_ws:
    :return:
    """
    try:
        # Use BL6:Chop:Skf1:SpeedUserReq to get chopper frequency to check only applied to chopper @ 30 Hz
        chopper_frequency = original_ws.run().getProperty('BL6:Chop:Skf1:SpeedUserReq').value.mean()
        if abs(chopper_frequency - 30) > 5.:
            # if chopper is not 30 Hz
            return output_ws

        # Check sample logs
        lead_max_wl = original_ws.run().getProperty('lead_max').value
        skip_min_wl = original_ws.run().getProperty('skip_min').value
    except RuntimeError as run_err:
        if original_ws.getInstrument().getFullName() != 'EQ-SANS':
            raise RuntimeError('Either lead_max or skip_min is not in workspace: {}'.format(run_err))
        else:
            # workspace other than EQ-SANS does not necessarily to have these sample logs
            return output_ws
    # END-TRY

    # Mask a range of X-values
    output_ws = MaskBins(output_ws, XMin=lead_max_wl, XMax=skip_min_wl)

    return output_ws
