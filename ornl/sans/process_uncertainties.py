from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import SetUncertainties
import numpy


def set_init_uncertainties(input_ws, output_ws=None):
    """
    Set the initial uncertainty of a MatrixWorkspace
    Mantid algorithm SetUncertainties will be called to make sure
    1: set the uncertainty to square root of intensity
    2: make sure all zero uncertainties will be set to 1

    In case of output workspace is None, the input workspace will be
    replaced by output workspace.

    :exception RuntimeError: output workspace (string) is empty

    :param input_ws: Input workspace
    :param output_ws: Output workspace (workspace name or instance) or None
    :return: reference to output workspace
    """
    # Set output workspace
    if output_ws is None:
        output_ws_name = str(input_ws)
    else:
        output_ws_name = str(output_ws).strip()
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

    return output_ws
