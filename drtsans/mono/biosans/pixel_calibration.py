from drtsans.pixel_calibration import load_and_apply_pixel_calibration as base_function

__all__ = ['load_and_apply_pixel_calibration', ]


def load_and_apply_pixel_calibration(input_workspace, output_workspace=None):
    r"""
    Update pixel positions and heights, as well as front and back tube widths, for a specific
    double panel detector array.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Input workspace to update
    output_workspace: str
        Optional name of the output workspace. if :py:obj:`None`, the name of ``input_workspace`` is used, thus
        calibrating the pixelsof the input workspace.
    """
    return base_function(input_workspace, output_workspace=output_workspace, components=['detector1', 'wing_detector'])
