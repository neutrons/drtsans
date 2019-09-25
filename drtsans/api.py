# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py
from drtsans.settings import unique_workspace_dundername as uwd
# https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Minus-v1.html
# https://docs.mantidproject.org/nightly/algorithms/RebinToWorkspace-v1.html
from mantid.simpleapi import CloneWorkspace, CreateSingleValuedWorkspace, Minus, mtd, RebinToWorkspace

__all__ = ['subtract_background']


def subtract_background(input_workspace, background, scale=1.0, scale_error=0.0,
                        output_workspace=None):
    r"""
    Subtract a prepared background from a prepared sample.

    Perform a rebin if sample and background have different binning.

    **Mantid algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`
    :ref:`CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>`
    :ref:`Minus <algm-Minus-v1>`
    :ref:`Multiply <algm-Multiply-v1>`
    :ref:`RebinToWorkspace <algm-RebinToWorkspace-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Sample workspace.
    background: str, ~mantid.api.MatrixWorkspace
        Background workspace.
    scale: float
        Rescale background intensities by this multiplicative factor before
        subtraction from the sample.
    scale_error: float
        Uncertainty in scale factor
    output_workspace: str
        Name of the sample corrected by the background. If :py:obj:`None`, then
        ``input_workspace`` will be overwritten.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    # get handle on input_workspace
    input_workspace = mtd[str(input_workspace)]

    workspaces_to_delete = []

    # make the background match the input_workspace binning if possible
    if input_workspace.isHistogramData() and background.isHistogramData():
        # rebin the background to match the input_workspace
        # this will do nothing if the x-axis is already the same
        # put it in the output_workspace to save space in memory
        background_rebinned = RebinToWorkspace(WorkspaceToRebin=background,
                                               WorkspaceToMatch=input_workspace,
                                               OutputWorkspace=uwd())
    else:
        # subtract will check that the x-axis is identical
        # otherwise the subtraction will fail
        background_rebinned = CloneWorkspace(InputWorkspace=background, OutputWorkspace=uwd())
    workspaces_to_delete.append(str(background_rebinned))

    # need to create a special object if the uncertainty in the scale was specified
    if scale_error != 0.:
        scale = CreateSingleValuedWorkspace(DataValue=scale, ErrorValue=scale_error,
                                            OutputWorkspace=uwd())
        workspaces_to_delete.append(str(scale))

    background_rebinned *= scale

    Minus(LHSWorkspace=input_workspace,
          RHSWorkspace=background_rebinned,
          OutputWorkspace=output_workspace)

    for name in workspaces_to_delete:
        if mtd.doesExist(name):
            mtd.remove(name)

    return mtd[output_workspace]
