from mantid.simpleapi import mtd, SetUncertainties
import numpy


def set_init_uncertainties(input_workspace, output_workspace=None):
    """
    Set the initial uncertainty of a :py:obj:`~mantid.api.MatrixWorkspace`

    Mantid algorithm :ref:`SetUncertainties <algm-SetUncertainties-v1>` will be called to make sure
    1: set the uncertainty to square root of intensity
    2: make sure all zero uncertainties will be set to 1

    In case of output workspace is py:obj:`None`, the input workspace will be
    replaced by output workspace.

    :exception RuntimeError: output workspace (string) is empty

    **Mantid algorithms used:**
    :ref:`SetUncertainties <algm-SetUncertainties-v1>`

    Parameters
    ----------
    input_workspace: ~mantid.api.MatrixWorkspace
        Input workspace
    output_workspace: str
        Output workspace (workspace name or instance) or py:obj:`None` for in-place operation

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

    return output_ws
